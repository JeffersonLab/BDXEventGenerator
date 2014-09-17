/*******************************************************************************
*									       *
* mcf_NTuIOFiles.h -- Utilities to manipulate files within the MCFIO Gen.      *
*        				Ntuple schema                          *
*									       *
*	P. Lebrun, September 1995.					       *
*									       *
*******************************************************************************/
void    mcf_ComputeNTuOffsets(nTuDDL *ddl); 
void    mcf_ComputeNTuLengths(nTuDDL *ddl); 
void    mcf_ComposeDoth(descrGenNtuple *dNtu, char *filename);
