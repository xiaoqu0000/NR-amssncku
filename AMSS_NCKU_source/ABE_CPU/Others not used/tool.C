#include <fstream>
#include <string>
// #include<
using namespace std;
/*void printss(int * a,int * b,int *c){
	int a1 = *a;
	int b1 = *b;
	int c1 = *c;
	printf("%d,%d,%d\n",1,2,3);
	printf("%d,%d,%d\n",a1,b1,c1);
}*/
int main()
{
	ifstream fin;
	ofstream fout;
	fin.open("tool_input.txt");
	fout.open("tool_output.txt");

	// ifstream fin1;
	// fin1.open("input1.txt");
	char buf[20];
	char buf1[20];

	while (fin >> buf)
	{
		// fin1>>buf1;
		// fout<<"if("<<buf<<") cudaFree("<<buf<<");\n";

		// cudaMalloc((void**)&(Mh_ #), matrix_size * sizeof(double));
		// fout<<"cudaMalloc((void**)&(Mh_"<<buf<<"), matrix_size * sizeof(double));"<<endl;

		// cudaMemcpy(Mh_ #, #, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
		// fout<<"cudaMemcpy(Mh_ "<<buf<<","<<buf<<", matrix_size * sizeof(double), cudaMemcpyHostToDevice);\n";

		// cudaMemcpy(#, Mh_ #, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		// fout<<"cudaMemcpy("<<buf<<", Mh_ "<<buf<<", matrix_size * sizeof(double), cudaMemcpyDeviceToHost);\n";

		// if(cg->[buf][i] != cg_gpu->[buf][i]){is_match = false; break;}
		fout << "delta = cg->fgfs[" << buf << "][i] - cg_gpu->fgfs[" << buf << "][i];" << endl;
		fout << "if(delta >1e-12 || delta < -1e-12){is_match = false; break;}" << endl;
	}
	/*int para = 167;
	for(int i = para;i<para+68;++i){
		fout<<"cg->fgfs["<<i<<"], ";
	}*/

	/*int array[3] = {0,1,2};
	int * p = array;
	printss(p++,p++,p++);*/
	return 0;
}
