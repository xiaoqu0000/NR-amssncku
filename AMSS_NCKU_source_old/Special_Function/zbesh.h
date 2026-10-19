
#ifndef ZBESH_H
#define ZBESH_H

#ifdef fortran1
#define f_zbesj zbesj
#endif
#ifdef fortran2
#define f_zbesj ZBESJ
#endif
#ifdef fortran3
#define f_zbesj zbesj_
#endif

extern "C"
{
	int f_zbesj(double &, double &, double &, int &,
				int &, double &, double &, int &, int &);
}
#endif /* ZBESH_H */
