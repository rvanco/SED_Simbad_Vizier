import cython

cdef extern from "c_astrom2.c":
	double to_jsky(float fmag, float mag)
	int search_vega_filter(char *system, char *filtername, float *Lambda, float *dlambda, float *fmag)

def c_to_jsky(fmag, mag):
	return to_jsky(fmag, mag)

def c_search_vega_filter(str_system, str_filtername) :
	byte_system    = bytes(str_system, 'UTF-8')
	cdef char* system = byte_system

	byte_filtername    = bytes(str_filtername, 'UTF-8')
	cdef char* filtername = byte_filtername

	cdef float init_Lambda = 0
	cdef float init_dlambda = 0
	cdef float init_fmag = 0

	cdef float* Lambda = &init_Lambda
	cdef float* dlambda = &init_dlambda
	cdef float* fmag = &init_fmag


	search_vega_filter(system, filtername, Lambda, dlambda, fmag)
	if Lambda[0]==0.0 and dlambda[0]==0.0 and fmag[0]==0.0 :
		print(f"The combination of the system {str_system} and the filter {str_filtername}, doesn't exist in the data.")
		return Lambda[0], dlambda[0], fmag[0]
	else :
		return Lambda[0], dlambda[0], fmag[0]
