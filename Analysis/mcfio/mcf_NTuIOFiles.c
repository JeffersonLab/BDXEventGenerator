/*******************************************************************************
*									       *
* mcf_NTuIOFiles.c -- Utilities to manipulate files within the MCFIO Gen.      *
*        				Ntuple schema                          *
*									       *
*	P. Lebrun, September 1995.					       *
*									       *
*******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/param.h>
#include <limits.h>
#include <time.h>
#include <rpc/types.h>
#include <sys/types.h>
#include <rpc/xdr.h>
#include <unistd.h>
#include <ctype.h>
#include "mcf_nTupleDescript.h"
#include "mcf_xdr.h"
#include "mcfio_Dict.h"
#include "mcfio_Direct.h"
#include "mcf_NTuIOFiles.h"
#include "mcf_NTuIOUtils.h"
#include "mcf_ntubld_db.h"
#ifndef False
#define False 0
#endif
#ifndef True
#define True 1
#endif

extern char *VarTypesNamesF77[N_VAR_TYPES];
extern char *VarTypesNamesC[N_VAR_TYPES];

extern struct line_title_c line_title_c_;
extern struct header_c header_c_;
extern struct variable_c variable_c_;

/*
** Ntuple identifiers list, initialized here and in mcfio_Util1
*/
nTuDDL **NTuDDLList = NULL;
int NumOfNTuples = 0;
bool_t McfNTuPleSaveDecoding = True;

static char *makeStructName(char *title, int orgStyle);
static size_t nDatVariable(varGenNtuple *varTmp);
static size_t sizeVariable(varGenNtuple *varTmp);


/*
** Compute the offsets by writing a simple program, stand alone, that uses
** the d/s
*/
void    mcf_ComputeNTuOffsets(nTuDDL *ddl)    
{
     char tmpName[128], *tc, *tc1, *nameCom;
     int i, j, l, fd, firstIndexed, nDat;
     char filenameInclude[128], filenameProgram[128], filenameExec[128];
     char filenameData[128], nameMaxIndex[32];
     char line[256];
     void **ptrBegVar;
     varGenNtuple *varTmp;
     descrGenNtuple *dNTu;
     FILE *Ffp;
     
     dNTu =ddl->descrNtu;  
     for (i=0; i< dNTu->numVariables; i++) { 
        varTmp = dNTu->variables[i];
        varTmp->offset = 0;
     }
        
     memset(tmpName, 0, 127);
     tc = tmpName;
     sprintf(tc, "tmp_%s_XXXXXX", ddl->descrNtu->title);
/* this is a kludge - we create a temporary file, close it, and use the name */
     fd = mkstemp(tmpName);
     if ( fd < 0 ) {
         fprintf(stderr, 
        " Can not compose a tempoary name in mcf_ComputeOffsets!");
         return;
     }
     tc1 = tc;
     close(fd);
     sprintf(filenameInclude, "%s.h", tc1);
     sprintf(filenameProgram, "%s.c", tc1);
     sprintf(filenameData, "%s.dat", tc1);
     strcpy(filenameExec, tc1);
     mcf_ComposeDoth(ddl->descrNtu, filenameInclude);
/*
** Compose a little moronic program that establishes the addresses of all 
** variables. There might be a better way, though.. However, this ought to be
** safe. 
*/
     Ffp =  fopen( filenameProgram, "w");
     fprintf(Ffp, "#include <stdio.h>\n");
     fprintf(Ffp, "#include \"%s\"\n",filenameInclude);
     if (dNTu->orgStyle == PARALLEL_ARRAY_NTU)
         fprintf(Ffp, "#define NUM_VAR %d\n", (dNTu->numVariables+3));
     else 
         fprintf(Ffp, "#define NUM_VAR %d\n",
                (dNTu->numVariables + 3 + dNTu->maxMultiplicity) );
         
     nameCom = makeStructName(dNTu->title, dNTu->orgStyle);
     
     fprintf(Ffp, "%s_struct tmpStruct; \n", nameCom);
     fprintf(Ffp, "main(int argc, char **argv)\n");
     fprintf(Ffp, "{\n");
     fprintf(Ffp, "    void *ptrBegVar[NUM_VAR];\n");
     fprintf(Ffp, "    FILE *Ffp;\n");
     fprintf(Ffp, "    int i;\n");
     fprintf(Ffp, "\n");
     fprintf(Ffp, "    ptrBegVar[0] = (void *) &tmpStruct.version[0];\n");
     fprintf(Ffp,
      "    ptrBegVar[1] = (void *) &tmpStruct.%s;\n",dNTu->nameIndex);
     for(i=0, firstIndexed=-1; i<dNTu->numVariables; i++) {
           if (dNTu->variables[i]->isFixedSize == False) {
                  firstIndexed = i; break;
           }       
     }
     dNTu->firstIndexed = firstIndexed;             
     if (dNTu->orgStyle == PARALLEL_ARRAY_NTU) {
         for(i=0; i<dNTu->numVariables; i++) {
           varTmp =  dNTu->variables[i];
           /*
           ** Assume that all the variables are properly 
           ** defined at this stage (e..g, coming from a valid DDL dbin file)
           ** and in order
           */
           tc = line;
           if ((varTmp->numDim == 0) && (varTmp->isFixedSize == True)) 
               sprintf(tc,
                "    ptrBegVar[%d] = (void *) &tmpStruct.%s%n",
                    (i+2), varTmp->name, &l);
           else           
               sprintf(tc,
                "    ptrBegVar[%d] = (void *) tmpStruct.%s%n",
                    (i+2), varTmp->name, &l);
           tc+=l;
           fprintf(Ffp, "%s;\n", line);    
         }
         fprintf(Ffp,
      "    ptrBegVar[%d] = (void *) tmpStruct.fence;\n",dNTu->numVariables+2);
     } else {
         for(i=0; i<dNTu->numVariables; i++) {
           varTmp =  dNTu->variables[i];
           tc = line;
           if (varTmp->isFixedSize == True) {
              if (varTmp->numDim == 0)
                sprintf(tc,
               "    ptrBegVar[%d] = (void *) &tmpStruct.%s%n",
                    (i+2), varTmp->name, &l);
               else 
                sprintf(tc,
               "    ptrBegVar[%d] = (void *) tmpStruct.%s%n",
                    (i+2), varTmp->name, &l);
           } else {
               if (varTmp->numDim == 0)         
                   sprintf(tc,
                   "    ptrBegVar[%d] = (void *) &tmpStruct.var[0].%s%n",
                    (i+2), varTmp->name, &l);
               else 
                   sprintf(tc,
                   "    ptrBegVar[%d] = (void *) tmpStruct.var[0].%s%n",
                    (i+2), varTmp->name, &l);
           }                                     
           fprintf(Ffp, "%s;\n", line); 
         }
         tc1 = dNTu->nameIndex;
         strcpy(nameMaxIndex, tc1);
         l = strlen(tc1); 
         if (l > 26) {
            strncpy(nameMaxIndex, tc1, 26);
            sprintf(&nameMaxIndex[26],"_max");
         } else
            sprintf(nameMaxIndex, "%s_max", tc1);
         fprintf(Ffp,"    for (i=0; i<%s; i++) \n", nameMaxIndex);
         tc = line;
         if (firstIndexed != -1) {
             varTmp = dNTu->variables[firstIndexed];
             sprintf(tc,
             "       ptrBegVar[i+%d] = (void *) &tmpStruct.var[i].%s%n",
              (2+dNTu->numVariables), varTmp->name, &l); tc+=l;
              if (varTmp->numDim > 0) for (j=0; j<varTmp->numDim; j++, tc+=l) 
                                             sprintf(tc, "[0]%n", &l);   
              fprintf(Ffp, "%s;\n", line);
         }      
         fprintf(Ffp,
      "    ptrBegVar[%d] = (void *) tmpStruct.fence;\n",
              dNTu->numVariables+2+dNTu->maxMultiplicity);
     }
     fprintf(Ffp, " ");
     fprintf(Ffp,"    Ffp = fopen(\"%s\",\"w\");\n",filenameData);
     fprintf(Ffp,"    fwrite((void *) ptrBegVar, sizeof(void *),\
(size_t) NUM_VAR, Ffp);\n");
     fprintf(Ffp,"    fclose(Ffp);\n");
     fprintf(Ffp,"}\n");
     fclose(Ffp);
     free(nameCom);
     /*
     ** Now compile, link and load this exec, read the result
     */
     sprintf(line,"rm -f %s", filenameExec);
     system(line);
#ifdef _HPUX_SOURCE
     sprintf(line,"cc -Aa -D_HPUX_SOURCE -o %s %s", 
                            filenameExec, filenameProgram);
#else
     sprintf(line,"cc -o %s %s", filenameExec, filenameProgram);
#endif
     system(line);
     sprintf(line,"./%s", filenameExec);
     system(line);
     if (dNTu->orgStyle == PARALLEL_ARRAY_NTU) nDat = dNTu->numVariables+3;
        else nDat = dNTu->numVariables+3+dNTu->maxMultiplicity;
     if (firstIndexed == -1) nDat = dNTu->numVariables+3;   
     ptrBegVar = (void **) malloc (sizeof(void *) * (nDat));
     Ffp = fopen(filenameData, "r");
     fread((void *) ptrBegVar, sizeof(void *), (size_t) nDat, Ffp);
     fclose(Ffp);
     /*
     ** remove garbage files..
     */
     remove(filenameData); remove(filenameProgram); remove(filenameExec);
     remove(filenameInclude);
     /*
     ** Convert these addresses to offsets
     */
     dNTu->multOffset =  ((long) ptrBegVar[1] - (long) ptrBegVar[0]);
     if (dNTu->orgStyle == PARALLEL_ARRAY_NTU) {
         dNTu->fenceOffset =
          ((long) ptrBegVar[dNTu->numVariables+2] - (long) ptrBegVar[0]); 
         for (i=0; i< dNTu->numVariables; i++) 
            dNTu->variables[i]->offset = 
                   ((long) ptrBegVar[i+2] - (long) ptrBegVar[0]);
     } else {
         for (i=0; i< dNTu->numVariables; i++) {
            varTmp = dNTu->variables[i]; 
            if (varTmp->isFixedSize)
                varTmp->offset = 
                    ((long) ptrBegVar[i+2] - (long) ptrBegVar[0]);
            else 
                varTmp->offset = 
                   ((long) ptrBegVar[i+2] - (long)ptrBegVar[firstIndexed+2]);
         }
         if (dNTu->subOffset != NULL) free(dNTu->subOffset);
         dNTu->subOffset =
            (long *) malloc(sizeof(long) * dNTu->maxMultiplicity);
         if (firstIndexed != -1) {    
            for (i=0; i<dNTu->maxMultiplicity; i++) 
                dNTu->subOffset[i] =          
               ((long) ptrBegVar[i+2+dNTu->numVariables]  - 
               (long) ptrBegVar[0]);
         }          
         dNTu->fenceOffset =
            ((long) ptrBegVar[dNTu->numVariables+2+dNTu->maxMultiplicity]
                     - (long) ptrBegVar[0]); 
     }   
     free(ptrBegVar);
}

/*
** Compute the lengths for the XDR Array statements. It is assumed that the 
** NTUple descriptor is sorted, no blank variables.
*/
void    mcf_ComputeNTuLengths(nTuDDL *ddl)    
{
     int i, j, lastTmp, sameType;
     size_t nDat, sizeItem;
     varGenNtuple *var1, *var2;
     descrGenNtuple *dNTu;
     
     dNTu =ddl->descrNtu;
     if (dNTu->firstIndexed != -1) lastTmp = dNTu->firstIndexed;
         else lastTmp = dNTu->numVariables;
     /*
     ** fixed size first.. 
     */   
     for (i=0; i<lastTmp; i++)
         dNTu->variables[i]->lengthW = nDatVariable(dNTu->variables[i]); 
/*
** This, in principle, is the optimized version, where we collaps single
** fields of the same type into an array. However, this is machine 
** dependant.
*/                                    
     for (i=0; i<lastTmp; i++) {
         var1 = dNTu->variables[i];
         if (var1->lengthW != 0) {
            nDat = nDatVariable(var1);
            j=i+1;
            sizeItem = sizeVariable(var1); 
            sameType = True;
            while ((j<lastTmp) && (sameType)) { 
               var2 = dNTu->variables[j];
               if (var2->type != var1->type) sameType = False;
               if (sameType  && ((( var2->offset - 
                                    var1->offset)/sizeItem) ==
                                     nDat)) {
                  nDat += nDatVariable(var2);
                  var2->lengthW = 0; j++;
               }
            }
            var1->lengthW = nDat;
            var1->lengthB = nDat*sizeItem;
         }
     } 
     /*
     ** The variable size, similar code. This fill is very simple if the 
     ** if the organisation is parallel arrays, as we can not implmenent
     ** compaction 
     */           
     if (dNTu->firstIndexed == -1) return;
     if (dNTu->orgStyle == PARALLEL_ARRAY_NTU) {
         for (i=dNTu->firstIndexed; i<dNTu->numVariables; i++) { 
              dNTu->variables[i]->lengthW
                          =  nDatVariable(dNTu->variables[i]);
              dNTu->variables[i]->lengthB = dNTu->variables[i]->lengthW   
                         * sizeVariable(dNTu->variables[i]);
         }                
     } else {
         for (i=dNTu->firstIndexed; i<dNTu->numVariables; i++) 
            dNTu->variables[i]->lengthW =  nDatVariable(dNTu->variables[i]);
         for (i=dNTu->firstIndexed; i<dNTu->numVariables; i++) {
             var1 = dNTu->variables[i];
             if (var1->lengthW != 0) {
                 nDat = nDatVariable(var1);
                 j=i+1;
                 sizeItem = sizeVariable(var1); 
                 sameType = True;
                 while ((j<dNTu->numVariables) && (sameType)) { 
                    var2 = dNTu->variables[j];
                    if (var2->type != var1->type) sameType = False;
                    if (sameType  && (((var2->offset - 
                                   var1->offset)/sizeItem) ==
                                     nDat)) {
                     nDat += nDatVariable(var2);
                     var2->lengthW = 0; j++; 
                  }
                }
                var1->lengthW = nDat;
                var1->lengthB = nDat*sizeItem;
             }
         }
     }
     
}
/*
** Compute, in size_t units (bytes, I hope) the length of a particular
** variable.  Only the fixed size part, we will have to multiplity 
** by the multiplicty in the XDR filter.
*/

static size_t nDatVariable(varGenNtuple *var)
{
   size_t n;
   int i;
    
    n=1;
    for (i=0; i<var->numDim; i++) n = n * var->dimensions[i];
    return n;
}
static size_t sizeVariable(varGenNtuple *var)
{
   size_t n;
   
   switch (var->type) {
        case BYTE_NTU: case CHARACTER_NTU:
           n = sizeof(char);
           break;
        case INTEGER2_NTU:
           n = sizeof(short);
           break;
        case LOGICAL_NTU: case INTEGER_NTU:
           n = sizeof(int);
           break;
        case REAL_NTU:
           n = sizeof(float);
           break;
        case DBL_PRECISION_NTU:
           n = sizeof(double);
           break;
        case COMPLEX_NTU:
           n = 2 * sizeof(float);
           break;
        case DBL_COMPLEX_NTU:
           n = 2 * sizeof(double);
           break;
        case POINTER_NTU:
           n = sizeof(void *);
           break;
        default : 
           fprintf(stderr, " mcf_ComputNTuLength, internal error \n");
           n = 0;
           break;
    }
    return n;
}

/*
** Compose the .h file. Called from NTuBldMenu and this file. The structure 
** is assumed valid. 
*/
void    mcf_ComposeDoth(descrGenNtuple *dNTu, char *filename)
{
    char *nameCom, line[FILENAME_MAX+500], *tmp, *version, *text, *tc, *tc2;
    char nameMaxIndex[32], nameTmpIndex[32];
    char nullDescr[4], *descrTmp;
    int i, j, l, kmode, nc, ncTot, nl, iv;
    time_t clock;
    FILE *Ffp;
    varGenNtuple *var;
    
    nameCom = makeStructName(dNTu->title, dNTu->orgStyle);
    strcpy(nullDescr, "? ");
    strcpy(line, filename);
    tc = strchr(line, '.');
    if (tc == NULL) {
         l = strlen(filename);
         tc = line; tc+=l;
    }     
    strcpy(tc,".h");   
    Ffp = fopen(line, "w");
    fprintf(Ffp,"/* ntuBuild\n");
    time(&clock);
    tmp = line; sprintf(tmp,"** Creation Date : %n", &l); tmp += l;
    strncpy(tmp,ctime(&clock), 24); tmp += 24; *tmp='\n'; tmp++; *tmp = '\0';
    fprintf(Ffp,line);
    fprintf(Ffp,"**  User Comments\n");
    text = dNTu->description;
    tc = text; 
    if (*tc == '\0') 
       fprintf(Ffp,"** no user comments\n");
    else {
       ncTot = strlen(tc); nc =0;
       while (nc < ncTot) {
            tc2 = strchr(tc,'\n');
            nl = (int) (tc2-tc)/sizeof(char);
            if ((tc2 == NULL) || (nl > 75)) nl = 75;
            strncpy(line, tc, nl); line[nl] = '\0';
            fprintf (Ffp,"** %s\n", line);
            tc += nl; nc += nl;
            if (*tc == '\n') {
               tc++;
               nc++;
            }
       }
    }
    fprintf(Ffp,"*/    \n");
    version = dNTu->version;
    text = dNTu->nameIndex;
    strcpy(nameTmpIndex, text);
    l = strlen(text); 
    if (l > 26) {
            strncpy(nameMaxIndex, text, 26);
            sprintf(&nameMaxIndex[26],"_max");
    } else
            sprintf(nameMaxIndex, "%s_max", text);
    fprintf(Ffp,"#define %s %d\n", nameMaxIndex, dNTu->maxMultiplicity);        
    if (dNTu->orgStyle == PARALLEL_ARRAY_NTU) {
        fprintf(Ffp, "typedef struct _%s_struct {\n", nameCom);
        /*
        ** The first 64 bits contain the version token, as a char[8] string
        ** floowed by the multiplicty variable, followed by an integer pad
        */ 
        fprintf(Ffp,"    char version[8]; /* Version token */\n");
        fprintf(Ffp,
        "    int %s; /* Generalized Ntuple Multiplicity value */ \n",
                                              nameTmpIndex);
        fprintf(Ffp,
        "    int padding; /* Padding for 64 bit architecture */ \n");
        for (iv=0; iv< dNTu->numVariables; iv++) {
           for (j=0; j<dNTu->numAvailable; j++)  
                if (dNTu->varOrdering[j] == iv)  i = j; 
           var = dNTu->variables[i];
           kmode = 0; if (var->isFixedSize != True) kmode = 1;
           if (var->description == NULL) descrTmp = nullDescr;
                  else descrTmp = var->description;
           tc = line;
           if ((var->type != COMPLEX_NTU) &&
               (var->type != DBL_COMPLEX_NTU)) { 
                  sprintf(tc,"    %s %n", VarTypesNamesC[var->type], &l);
                  tc +=l;            
                  if ((var->numDim == 0) && (kmode ==0))
                            sprintf(tc," %s; /* %s */",
                                              var->name, descrTmp);
                  else if (var->numDim == 0) {
                        sprintf(tc," %s[%s]; /* %s */",
                        var->name, nameMaxIndex, descrTmp); 
                  } else { 
                       sprintf(tc," %s%n",var->name, &l); tc+=l;
                       if (kmode == 1) {
                           sprintf(tc, "[%s]%n", nameMaxIndex, &l);
                           tc +=l;
                   }    
                   for (j=var->numDim-1; j>-1; j--, tc+=l) 
                             sprintf(tc,"[%d]%n", var->dimensions[j], &l); 
                     
                        sprintf (tc,"; /* %s */", descrTmp);
                  }  
               } else { /* got to convert to float or dbl */
                    if (var->type == COMPLEX_NTU)
                       sprintf(tc,"    float %n", &l);
                       
                     else if (var->type == DBL_COMPLEX_NTU)
                       sprintf(tc,"    double %n", &l);
                       
                     tc +=l;            
                     if ((var->numDim == 0) && (kmode ==0)) 
                          sprintf(tc," %s[2]; /* %s */", var->name, descrTmp);
                     else if (var->numDim == 0) {
                       sprintf(tc," %s[%s][2]; /* %s */",
                        var->name, nameMaxIndex, descrTmp); 
                     } else { 
                       sprintf(tc," %s%n",var->name, &l); tc+=l;
                       if (kmode == 1) {
                          sprintf(tc, "[%s]%n", nameMaxIndex, &l);
                         tc +=l;
                       }    
                       for (j=var->numDim-1; j>-1; j--, tc+=l) 
                             sprintf(tc,"[%d]%n", var->dimensions[j], &l); 
                       sprintf (tc,"[2]; /* %s */", descrTmp);
                    }
               }
               fprintf(Ffp,"%s\n", line);
           }
        fprintf(Ffp,"    int fence[2]; \n");
        fprintf(Ffp,"} %s_struct; \n", nameCom);
     }else { 
     /*
     ** The other type of organisation, using structure
     */
        fprintf(Ffp, "typedef struct _%s_v_struct{\n", nameCom);
        for (iv=0; iv< dNTu->numVariables; iv++) {
           for (j=0; j<dNTu->numAvailable; j++)  
                if (dNTu->varOrdering[j] == iv)  i = j; 
           var = dNTu->variables[i];
           if (var->isFixedSize == False) {
               tc = line;
               if (var->type == COMPLEX_NTU)
                       sprintf(tc,"    float %n", &l);
               else if (var->type == DBL_COMPLEX_NTU)
                       sprintf(tc,"    double %n", &l);
               else 
                     sprintf(tc,"    %s %n", VarTypesNamesC[var->type], &l);
               tc +=l;            
               sprintf(tc," %s%n",var->name, &l); tc+=l;
               if (var->numDim != 0) {
                    for (j=var->numDim-1; j>-1; j--, tc+=l) 
                           sprintf(tc,"[%d]%n", var->dimensions[j], &l); 
               }
               if ((var->type == COMPLEX_NTU) ||
                   (var->type == DBL_COMPLEX_NTU)) {
                       sprintf (tc,"[2]%n",&l);
                       tc += l;
               }       
               if (var->description == NULL) descrTmp = nullDescr;
                    else descrTmp = var->description;
               sprintf(tc,"; /* %s */%n", descrTmp, &l); tc += l;
               fprintf(Ffp,"%s\n", line);
            }
        }
        fprintf(Ffp,"} %s_v_struct; \n", nameCom);
        fprintf(Ffp,"/* ----- */  \n");
        /*
        ** the mother structure now
        */
        fprintf(Ffp, "typedef struct _%s_struct{\n", nameCom);
        fprintf(Ffp,"    char version[8]; /* Version token */\n");
        fprintf(Ffp,
        "    int %s; /* Generalized Ntuple Multiplicity value */ \n",
                                              nameTmpIndex);
        fprintf(Ffp,
        "    int padding; /* Padding for 64 bit architecture */ \n");
        for (iv=0; iv< dNTu->numVariables; iv++) {
           for (j=0; j<dNTu->numAvailable; j++)  
                if (dNTu->varOrdering[j] == iv)  i = j; 
           var = dNTu->variables[i];
           if (var->isFixedSize == True) {
              tc = line;
               if (var->type == COMPLEX_NTU)
                       sprintf(tc,"    float %n", &l);
               else if (var->type == DBL_COMPLEX_NTU)
                       sprintf(tc,"    double %n", &l);
               else 
                     sprintf(tc,"    %s %n", VarTypesNamesC[var->type], &l);
               tc +=l;            
               sprintf(tc," %s%n",var->name, &l); tc+=l;
               if (var->numDim != 0) {
                    for (j=var->numDim-1; j>-1; j--, tc+=l) 
                           sprintf(tc,"[%d]%n", var->dimensions[j], &l); 
               }
               if ((var->type == COMPLEX_NTU) ||
                     (var->type == DBL_COMPLEX_NTU)) {
                       sprintf (tc,"[2]%n",&l);
                       tc += l;
               }       
               if (var->description == NULL) descrTmp = nullDescr;
                    else descrTmp = var->description;
               sprintf(tc,"; /* %s */%n", descrTmp, &l); tc += l;
               fprintf(Ffp,"%s\n", line);
            }
        }
        fprintf(Ffp,
        "    %s_v_struct var[%s]; /* The array of substructures */\n",  
                              nameCom, nameMaxIndex); 
        fprintf(Ffp,"    int fence[2]; \n");
        fprintf(Ffp,"} %s_struct; \n", nameCom);
    }    
    free(nameCom);
    fclose(Ffp);
    
}

void mcfioC_SetForSaveDecoding(int val)
{
     if(val != 0) McfNTuPleSaveDecoding = True;
     else McfNTuPleSaveDecoding = False;
}
     
static char *makeStructName(char *title, int orgStyle)
{
    char *out;
    int i, l, nMax;
    
    l = strlen(title);
    if (orgStyle  == PARALLEL_ARRAY_NTU) nMax = 23;  
    else nMax = 21;
    if (l > nMax) l = nMax;
    out = (char *) malloc(sizeof(char) * (l+1));
    strncpy(out, title, l); out[l]='\0';    
    for (i=0; i<l; i++) if (out[i] == ' ') out[i] = '_';
    return out;
} 
