
#ifdef newc
#include <cstdio>
using namespace std;
#else
#include <stdio.h>
#endif

#include "checkpoint.h"
#include "misc.h"
#include "fmisc.h"
#include "parameters.h"
#include "data_motion_helper.h"

namespace {

constexpr int kCheckpointNameLength = 100;
constexpr int kParameterLineLength = 256;

void abort_open_failure(const char* file_name, const char* action) {
	cout << "Can't open " << file_name << " for check point " << action << "." << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
}

void open_checkpoint_output(ofstream& outfile, const char* file_name) {
	outfile.open(file_name, ios::out | ios::trunc);
	if (!outfile) {
		abort_open_failure(file_name, "out");
	}
}

void open_checkpoint_input(ifstream& infile, const char* file_name) {
	infile.open(file_name);
	if (!infile) {
		abort_open_failure(file_name, "in");
	}
}

void append_variable_node(MyList<var>*& check_list, var* variable) {
	if (!variable) {
		return;
	}

	if (check_list) {
		check_list->insert(variable);
	} else {
		check_list = new MyList<var>(variable);
	}
}

string resolve_output_dir() {
	map<string, string>::iterator output_dir_it = parameters::str_par.find("output dir");
	if (output_dir_it != parameters::str_par.end()) {
		return output_dir_it->second;
	}

	char parameter_line[kParameterLineLength];
	string line_text;
	string group_name;
	string key_name;
	string value_text;
	int section_index = 0;
	char parameter_file[50];

	map<string, string>::iterator input_file_it = parameters::str_par.find("inputpar");
	if (input_file_it != parameters::str_par.end()) {
		strcpy(parameter_file, (input_file_it->second).c_str());
	} else {
		cout << "Error inputpar" << endl;
		exit(0);
	}

	ifstream input_stream(parameter_file, ifstream::in);
	if (!input_stream.good()) {
		cout << "Can not open parameter file " << parameter_file << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	string output_dir;
	for (int line_number = 1; input_stream.good(); ++line_number) {
		input_stream.getline(parameter_line, kParameterLineLength);
		line_text = parameter_line;

		const int parse_status = misc::parse_parts(line_text, group_name, key_name, value_text, section_index);
		if (parse_status == -1) {
			cout << "error reading parameter file " << parameter_file << " in line " << line_number << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (parse_status == 0) {
			continue;
		}

		if (group_name == "ABE" && key_name == "output dir") {
			output_dir = value_text;
		}
	}
	input_stream.close();
	parameters::str_par.insert(map<string, string>::value_type("output dir", output_dir));
	return output_dir;
}

} 

checkpoint::checkpoint(bool checked, const char fname[], int myrank) : filename(0), CheckList(0), checkedrun(checked)
{
	out_dir = resolve_output_dir();

	I_Print = (myrank == 0);

	int i = strlen(fname);
	filename = new char[i + 30];

#ifdef CHECKDETAIL
	char cmd[80];
	if (!checkedrun) {
		sprintf(cmd, "rm -rf %s/%d", out_dir.c_str(), myrank);
		system(cmd);
		sprintf(cmd, "mkdir %s/%d", out_dir.c_str(), myrank);
		system(cmd);
	}
	sprintf(filename, "%s/%d/%s", out_dir.c_str(), myrank, fname);
#else
	sprintf(filename, "%s/%s", out_dir.c_str(), fname);
#endif
	if (myrank == 0) {
		cout << " checkpoint class created " << endl;
	}
}
checkpoint::~checkpoint()
{
	CheckList->clearList();
	if (I_Print)
		delete[] filename;
}

void checkpoint::addvariable(var *VV)
{
	append_variable_node(CheckList, VV);
}
void checkpoint::addvariablelist(MyList<var> *VL)
{
	while (VL)
	{
		append_variable_node(CheckList, VL->data);
		VL = VL->next;
	}
}
#ifndef CHECKDETAIL
void checkpoint::writecheck_cgh(double time, cgh *GH)
{
	ofstream outfile;

	if (I_Print)
	{
		char fname[kCheckpointNameLength];
		sprintf(fname, "%s_cgh.CHK", filename);

		open_checkpoint_output(outfile, fname);
		outfile.write((char *)&time, sizeof(double));
		outfile.write((char *)&(GH->levels), sizeof(int));
		outfile.write((char *)&(GH->movls), sizeof(int));
		outfile.write((char *)&(GH->BH_num_in), sizeof(int));
		outfile.write((char *)GH->grids, GH->levels * sizeof(int));
		outfile.write((char *)GH->Lt, GH->levels * sizeof(double));
		for (int lev = 0; lev < GH->levels; lev++)
		{
			for (int grd = 0; grd < GH->grids[lev]; grd++)
			{
				outfile.write((char *)GH->bbox[lev][grd], 6 * sizeof(double));
				outfile.write((char *)GH->shape[lev][grd], 3 * sizeof(int));
				outfile.write((char *)GH->handle[lev][grd], 3 * sizeof(double));
			}
			for (int ibh = 0; ibh < GH->BH_num_in; ibh++)
			{
				outfile.write((char *)GH->Porgls[lev][ibh], 3 * sizeof(double));
			}
		}
	}
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	// write variable data
	for (int lev = 0; lev < GH->levels; lev++)
	{
		MyList<Patch> *PL = GH->PatL[lev];
		DataMotionHelper::pull_all_data_to_host(PL, myrank, CheckList);
		while (PL)
		{
			Patch *PP = PL->data;
			int nn = PP->shape[0] * PP->shape[1] * PP->shape[2];
			MyList<var> *VL = CheckList;
			while (VL)
			{
				double *databuffer = Parallel::Collect_Data(PP, VL->data);
				if (I_Print)
					outfile.write((char *)databuffer, sizeof(double) * nn);
				if (databuffer)
					delete[] databuffer;
				VL = VL->next;
			}
			PL = PL->next;
		}
	}

	if (I_Print)
		outfile.close();
}
void checkpoint::readcheck_cgh(double &time, cgh *GH, int myrank, int nprocs, int Symmetry)
{
	int DIM = dim;
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_cgh.CHK", filename);

	open_checkpoint_input(infile, fname);

	infile.seekg(0, ios::beg);
	infile.read((char *)&time, sizeof(double));
	if (I_Print)
		cout << "check cgh in at t = " << time << endl;
	infile.read((char *)&(GH->levels), sizeof(int));
	infile.read((char *)&(GH->movls), sizeof(int));
	infile.read((char *)&(GH->BH_num_in), sizeof(int));
	GH->grids = new int[GH->levels];
	GH->bbox = new double **[GH->levels];
	GH->shape = new int **[GH->levels];
	GH->handle = new double **[GH->levels];
	GH->PatL = new MyList<Patch> *[GH->levels];
	GH->Lt = new double[GH->levels];
	GH->Porgls = new double **[GH->levels];
#if (RPB == 1)
	GH->bdsul = new MyList<Parallel::pointstru_bam> *[GH->levels];
	GH->rsul = new MyList<Parallel::pointstru_bam> *[GH->levels];
#endif
	infile.read((char *)GH->grids, GH->levels * sizeof(int));
	infile.read((char *)GH->Lt, GH->levels * sizeof(double));
	for (int lev = 0; lev < GH->levels; lev++)
	{
		GH->bbox[lev] = new double *[GH->grids[lev]];
		GH->shape[lev] = new int *[GH->grids[lev]];
		GH->handle[lev] = new double *[GH->grids[lev]];
		GH->Porgls[lev] = new double *[GH->BH_num_in];
		for (int grd = 0; grd < GH->grids[lev]; grd++)
		{
			GH->bbox[lev][grd] = new double[6];
			GH->shape[lev][grd] = new int[3];
			GH->handle[lev][grd] = new double[3];
			infile.read((char *)GH->bbox[lev][grd], 6 * sizeof(double));
			infile.read((char *)GH->shape[lev][grd], 3 * sizeof(int));
			infile.read((char *)GH->handle[lev][grd], 3 * sizeof(double));
		}
		for (int ibh = 0; ibh < GH->BH_num_in; ibh++)
		{
			GH->Porgls[lev][ibh] = new double[dim];
			infile.read((char *)GH->Porgls[lev][ibh], 3 * sizeof(double));
		}
	}

	for (int lev = 0; lev < GH->levels; lev++)
		GH->PatL[lev] = GH->construct_patchlist(lev, Symmetry);

	GH->compose_cgh(nprocs);
	// write variable data
	for (int lev = 0; lev < GH->levels; lev++)
	{
		MyList<Patch> *PL = GH->PatL[lev];
		while (PL)
		{
			Patch *PP = PL->data;
			int nn = PP->shape[0] * PP->shape[1] * PP->shape[2];
			double *databuffer = new double[nn];
			MyList<var> *VL = CheckList;
			while (VL)
			{
				infile.read((char *)databuffer, sizeof(double) * nn);

				{
					MyList<Block> *BL = PP->blb;
					while (BL)
					{
						Block *cg = BL->data;
						if (myrank == cg->rank)
						{
							f_copy(DIM, cg->bbox, cg->bbox + DIM, cg->shape, cg->fgfs[VL->data->sgfn],
								   PP->bbox, PP->bbox + DIM, PP->shape, databuffer,
								   cg->bbox, cg->bbox + DIM);
						}
						if (BL == PP->ble)
							break;
						BL = BL->next;
					}
				}

				VL = VL->next;
			}
			delete[] databuffer;
			PL = PL->next;
		}
	}

	infile.close();
}
void checkpoint::writecheck_sh(double time, ShellPatch *SH)
{
	ofstream outfile;

	if (I_Print)
	{
		char fname[kCheckpointNameLength];
		sprintf(fname, "%s_sh.CHK", filename);

		open_checkpoint_output(outfile, fname);
		outfile.write((char *)&time, sizeof(double));
	}

	// write variable data
	MyList<ss_patch> *Pp = SH->PatL;
	while (Pp)
	{
		int nn = Pp->data->shape[0] * Pp->data->shape[1] * Pp->data->shape[2];
		MyList<var> *VL = CheckList;
		while (VL)
		{
			double *databuffer = SH->Collect_Data(Pp->data, VL->data);
			if (I_Print)
				outfile.write((char *)databuffer, sizeof(double) * nn);
			if (databuffer)
				delete[] databuffer;
			VL = VL->next;
		}
		Pp = Pp->next;
	}

	if (I_Print)
		outfile.close();
}
void checkpoint::readcheck_sh(ShellPatch *SH, int myrank)
{
	int DIM = dim;
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_sh.CHK", filename);

	open_checkpoint_input(infile, fname);

	double time;
	infile.seekg(0, ios::beg);
	infile.read((char *)&time, sizeof(double));
	if (I_Print)
		cout << "check ShellPatch in at t = " << time << endl;

	// because we assume the shell patch is fixed we can leave the composing to other routine

	MyList<ss_patch> *Pp = SH->PatL;
	while (Pp)
	{
		int nn = Pp->data->shape[0] * Pp->data->shape[1] * Pp->data->shape[2];
		double *databuffer = new double[nn];
		MyList<var> *VL = CheckList;
		while (VL)
		{
			infile.read((char *)databuffer, sizeof(double) * nn);
			MyList<Block> *BL = Pp->data->blb;
			while (BL)
			{
				Block *cg = BL->data;
				if (myrank == cg->rank)
				{
					f_copy(DIM, cg->bbox, cg->bbox + DIM, cg->shape, cg->fgfs[VL->data->sgfn],
						   Pp->data->bbox, Pp->data->bbox + DIM, Pp->data->shape, databuffer,
						   cg->bbox, cg->bbox + DIM);
				}
				if (BL == Pp->data->ble)
					break;
				BL = BL->next;
			}
			VL = VL->next;
		}
		delete[] databuffer;
		Pp = Pp->next;
	}

	infile.close();
}
void checkpoint::write_Black_Hole_position(int BH_num_input, int BH_num, double **Porg0, double **Porgbr, double *Mass)
{
	ofstream outfile;

	if (I_Print)
	{
		char fname[kCheckpointNameLength];
		sprintf(fname, "%s_BHp.CHK", filename);

		open_checkpoint_output(outfile, fname);
		outfile.write((char *)&BH_num_input, sizeof(int));
		outfile.write((char *)&BH_num, sizeof(int));
		outfile.write((char *)Mass, 3 * sizeof(double));
		for (int i = 0; i < BH_num; i++)
		{
			outfile.write((char *)Porg0[i], 3 * sizeof(double));
			outfile.write((char *)Porgbr[i], 3 * sizeof(double));
		}

		outfile.close();
	}
}
void checkpoint::read_Black_Hole_position(int &BH_num_input, int &BH_num, double **&Porg0, double *&Pmom,
										  double *&Spin, double *&Mass, double **&Porgbr, double **&Porg,
										  double **&Porg1, double **&Porg_rhs)
{
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_BHp.CHK", filename);

	open_checkpoint_input(infile, fname);
	if (I_Print)
		cout << "checking in Black_Hole_position" << endl;

	infile.seekg(0, ios::beg);
	infile.read((char *)&BH_num_input, sizeof(int));
	infile.read((char *)&BH_num, sizeof(int));
	// these arrays will be deleted when bssn_class is deleted
	Pmom = new double[3 * BH_num];
	Spin = new double[3 * BH_num];
	Mass = new double[BH_num];
	Porg0 = new double *[BH_num];
	Porgbr = new double *[BH_num];
	Porg = new double *[BH_num];
	Porg1 = new double *[BH_num];
	Porg_rhs = new double *[BH_num];
	infile.read((char *)Mass, 3 * sizeof(double));
	for (int i = 0; i < BH_num; i++)
	{
		Porg0[i] = new double[3];
		Porgbr[i] = new double[3];
		Porg[i] = new double[3];
		Porg1[i] = new double[3];
		Porg_rhs[i] = new double[3];
		infile.read((char *)Porg0[i], 3 * sizeof(double));
		infile.read((char *)Porgbr[i], 3 * sizeof(double));
	}

	infile.close();
}
void checkpoint::write_bssn(double LastDump, double Last2dDump, double LastAnas)
{
	ofstream outfile;

	if (I_Print)
	{
		char fname[kCheckpointNameLength];
		sprintf(fname, "%s_bssn.CHK", filename);

		open_checkpoint_output(outfile, fname);
		outfile.write((char *)&LastDump, sizeof(double));
		outfile.write((char *)&Last2dDump, sizeof(double));
		outfile.write((char *)&LastAnas, sizeof(double));

		outfile.close();
	}
}
void checkpoint::read_bssn(double &LastDump, double &Last2dDump, double &LastAnas)
{
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_bssn.CHK", filename);

	open_checkpoint_input(infile, fname);
	if (I_Print)
		cout << "checking in bssn parameters" << endl;

	infile.seekg(0, ios::beg);
	infile.read((char *)&LastDump, sizeof(double));
	infile.read((char *)&Last2dDump, sizeof(double));
	infile.read((char *)&LastAnas, sizeof(double));

	infile.close();
}
#else
void checkpoint::write_bssn(double LastDump, double Last2dDump, double LastAnas)
{
	ofstream outfile;

	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_bssn.CHK", filename);

	open_checkpoint_output(outfile, fname);
	outfile.setf(ios::scientific, ios::floatfield);
	outfile.precision(16);
	outfile << LastDump << " ";
	outfile << Last2dDump << " ";
	outfile << LastAnas << " " << endl;

	outfile.close();
}
void checkpoint::read_bssn(double &LastDump, double &Last2dDump, double &LastAnas)
{
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_bssn.CHK", filename);

	open_checkpoint_input(infile, fname);
	if (I_Print)
		cout << "checking in bssn parameters" << endl;

	infile.seekg(0, ios::beg);
	infile >> LastDump;
	infile >> Last2dDump;
	infile >> LastAnas;

	infile.close();
}
void checkpoint::write_Black_Hole_position(int BH_num_input, int BH_num, double **Porg0, double **Porgbr)
{
	ofstream outfile;

	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_BHp.CHK", filename);

	open_checkpoint_output(outfile, fname);
	outfile.setf(ios::scientific, ios::floatfield);
	outfile.precision(16);
	outfile << BH_num_input << " ";
	outfile << BH_num << " ";
	for (int i = 0; i < BH_num; i++)
	{
		for (int j = 0; j < 3; j++)
			outfile << Porg0[i][j] << " ";
		for (int j = 0; j < 3; j++)
			outfile << Porgbr[i][j] << " ";
	}

	outfile << endl;
	outfile.close();
}
void checkpoint::read_Black_Hole_position(int &BH_num_input, int &BH_num, double **&Porg0, double *&Pmom,
										  double *&Spin, double *&Mass, double **&Porgbr, double **&Porg,
										  double **&Porg1, double **&Porg_rhs)
{
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_BHp.CHK", filename);

	open_checkpoint_input(infile, fname);
	if (I_Print)
		cout << "checking in Black_Hole_position" << endl;

	infile.seekg(0, ios::beg);
	infile >> BH_num_input;
	infile >> BH_num;
	// these arrays will be deleted when bssn_class is deleted
	Pmom = new double[3 * BH_num];
	Spin = new double[3 * BH_num];
	Mass = new double[BH_num];
	Porg0 = new double *[BH_num];
	Porgbr = new double *[BH_num];
	Porg = new double *[BH_num];
	Porg1 = new double *[BH_num];
	Porg_rhs = new double *[BH_num];
	for (int i = 0; i < BH_num; i++)
	{
		Porg0[i] = new double[3];
		Porgbr[i] = new double[3];
		Porg[i] = new double[3];
		Porg1[i] = new double[3];
		Porg_rhs[i] = new double[3];
		for (int j = 0; j < 3; j++)
			infile >> Porg0[i][j];
		for (int j = 0; j < 3; j++)
			infile >> Porgbr[i][j];
	}

	infile.close();
}
void checkpoint::writecheck_cgh(double time, cgh *GH)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	ofstream outfile;

	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_cgh.CHK", filename);

	open_checkpoint_output(outfile, fname);
	outfile.setf(ios::scientific, ios::floatfield);
	outfile.precision(16);
	outfile << time << " ";
	outfile << (GH->levels) << " ";
	outfile << (GH->movls) << " ";
	outfile << (GH->BH_num_in) << " ";
	for (int j = 0; j < GH->levels; j++)
		outfile << GH->grids[j] << " ";
	for (int j = 0; j < GH->levels; j++)
		outfile << GH->Lt[j] << " ";
	for (int lev = 0; lev < GH->levels; lev++)
	{
		for (int grd = 0; grd < GH->grids[lev]; grd++)
		{
			for (int j = 0; j < 6; j++)
				outfile << GH->bbox[lev][grd][j] << " ";
			for (int j = 0; j < 3; j++)
				outfile << GH->shape[lev][grd][j] << " ";
			for (int j = 0; j < 3; j++)
				outfile << GH->handle[lev][grd][j] << " ";
		}
		for (int ibh = 0; ibh < GH->BH_num_in; ibh++)
		{
			for (int j = 0; j < 3; j++)
				outfile << GH->Porgls[lev][ibh][j] << " ";
		}
	}
	// write variable data
	for (int lev = 0; lev < GH->levels; lev++)
	{
		MyList<Patch> *PL = GH->PatL[lev];
		int cnt = 0;
		while (PL)
		{
			cnt++;
			PL = PL->next;
		}
		outfile << cnt << " ";
		PL = GH->PatL[lev];
		while (PL)
		{
			Patch *PP = PL->data;
			outfile << PP->lev << " ";
			for (int j = 0; j < 3; j++)
				outfile << PP->shape[j] << " ";
			for (int j = 0; j < 6; j++)
				outfile << PP->bbox[j] << " ";
			for (int j = 0; j < 3; j++)
				outfile << PP->lli[j] << " ";
			for (int j = 0; j < 3; j++)
				outfile << PP->uui[j] << " ";

			MyList<Block> *BP = PP->blb;
			cnt = 0;
			while (BP)
			{
				Block *cg = BP->data;
				cnt++;
				if (BP == PP->ble)
					break;
				BP = BP->next;
			}
			outfile << cnt << " ";

			BP = PP->blb;
			while (BP)
			{
				Block *cg = BP->data;
				for (int j = 0; j < 3; j++)
					outfile << cg->shape[j] << " ";
				for (int j = 0; j < 6; j++)
					outfile << cg->bbox[j] << " ";
				outfile << cg->rank << " " << cg->lev << " " << cg->cgpu << " "
						<< cg->ingfs << " " << cg->fngfs << " ";
				if (myrank == cg->rank)
				{
					MyList<var> *VL = CheckList;
					int NN = cg->shape[0] * cg->shape[1] * cg->shape[2];
					while (VL)
					{
						for (int j = 0; j < NN; j++)
							outfile << cg->fgfs[VL->data->sgfn][j] << " ";
						VL = VL->next;
					}
				}
				if (BP == PP->ble)
					break;
				BP = BP->next;
			}
			PL = PL->next;
		}
	}

	outfile << endl;
	outfile.close();
}
void checkpoint::readcheck_cgh(double &time, cgh *GH, int myrank, int nprocs, int Symmetry)
{
	int DIM = dim;
	ifstream infile;
	char fname[kCheckpointNameLength];
	sprintf(fname, "%s_cgh.CHK", filename);

	open_checkpoint_input(infile, fname);

	infile.seekg(0, ios::beg);
	infile >> time;
	if (I_Print)
		cout << "check cgh in at t = " << time << endl;
	infile >> (GH->levels);
	infile >> (GH->movls);
	infile >> (GH->BH_num_in);
	GH->grids = new int[GH->levels];
	GH->bbox = new double **[GH->levels];
	GH->shape = new int **[GH->levels];
	GH->handle = new double **[GH->levels];
	GH->PatL = new MyList<Patch> *[GH->levels];
	GH->Lt = new double[GH->levels];
	GH->Porgls = new double **[GH->levels];
#if (RPB == 1)
	GH->bdsul = new MyList<Parallel::pointstru_bam> *[GH->levels];
	GH->rsul = new MyList<Parallel::pointstru_bam> *[GH->levels];
#endif
	for (int j = 0; j < GH->levels; j++)
		infile >> GH->grids[j];
	for (int j = 0; j < GH->levels; j++)
		infile >> GH->Lt[j];
	for (int lev = 0; lev < GH->levels; lev++)
	{
		GH->bbox[lev] = new double *[GH->grids[lev]];
		GH->shape[lev] = new int *[GH->grids[lev]];
		GH->handle[lev] = new double *[GH->grids[lev]];
		GH->Porgls[lev] = new double *[GH->BH_num_in];
		for (int grd = 0; grd < GH->grids[lev]; grd++)
		{
			GH->bbox[lev][grd] = new double[6];
			GH->shape[lev][grd] = new int[3];
			GH->handle[lev][grd] = new double[3];
			for (int j = 0; j < 6; j++)
				infile >> GH->bbox[lev][grd][j];
			for (int j = 0; j < 3; j++)
				infile >> GH->shape[lev][grd][j];
			for (int j = 0; j < 3; j++)
				infile >> GH->handle[lev][grd][j];
		}
		for (int ibh = 0; ibh < GH->BH_num_in; ibh++)
		{
			GH->Porgls[lev][ibh] = new double[dim];
			for (int j = 0; j < 3; j++)
				infile >> GH->Porgls[lev][ibh][j];
		}
	}

	// read variable data
	for (int lev = 0; lev < GH->levels; lev++)
	{
		int cnt;
		infile >> cnt;
		GH->PatL[lev] = 0;

		MyList<Patch> *gp;
		// loop of patach
		for (int cj = 0; cj < cnt; cj++)
		{
			if (GH->PatL[lev])
			{
				gp->next = new MyList<Patch>;
				gp = gp->next;
			}
			else
			{
				GH->PatL[lev] = gp = new MyList<Patch>;
			}
			gp->data = new Patch();
			infile >> gp->data->lev;
			for (int j = 0; j < 3; j++)
				infile >> gp->data->shape[j];
			for (int j = 0; j < 6; j++)
				infile >> gp->data->bbox[j];
			for (int j = 0; j < 3; j++)
				infile >> gp->data->lli[j];
			for (int j = 0; j < 3; j++)
				infile >> gp->data->uui[j];
			gp->next = 0;
			gp->data->blb = 0;
			gp->data->ble = 0;
			// loop of Block
			int bnt;
			infile >> bnt;

			MyList<Block> *cg;
			for (int bj = 0; bj < bnt; bj++)
			{
				if (gp->data->blb)
				{
					cg->next = new MyList<Block>;
					cg = cg->next;
				}
				else
				{
					gp->data->blb = cg = new MyList<Block>;
				}
				double tbbox[6];
				int tshape[3];
				int trank, tlev, tcgpu, tingfs, tfngfs;
				for (int j = 0; j < 3; j++)
					infile >> tshape[j];
				for (int j = 0; j < 6; j++)
					infile >> tbbox[j];
				infile >> trank >> tlev >> tcgpu >> tingfs >> tfngfs;
				cg->data = new Block(dim, tshape, tbbox, trank, tingfs, tfngfs, tlev, tcgpu);
				cg->next = 0;
// if read fake check data, comment out this part
#if 1
				if (myrank == cg->data->rank)
				{
					MyList<var> *VL = CheckList;
					int NN = cg->data->shape[0] * cg->data->shape[1] * cg->data->shape[2];
					while (VL)
					{
						for (int j = 0; j < NN; j++)
							infile >> cg->data->fgfs[VL->data->sgfn][j];
						VL = VL->next;
					}
				}
#endif
			}
			gp->data->ble = cg;
		}

#if (RPB == 1)
		// we need distributed box of PatL[lev] and PatL[lev-1]
		if (lev > 0)
		{
			Parallel::Constr_pointstr_OutBdLow2Hi(PatL[lev], PatL[lev - 1], bdsul[lev]);
			Parallel::Constr_pointstr_Restrict(PatL[lev], PatL[lev - 1], rsul[lev]);
		}
		else
		{
			bdsul[lev] = 0;
			rsul[lev] = 0;
		}
#endif
	}

	infile.close();
}
#endif
