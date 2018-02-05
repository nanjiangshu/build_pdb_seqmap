/*
 * =====================================================================================
 *       Filename:  getseqresseq.cpp
 *
 *    Description:  Obtain the amino acid sequences from the SEQRES records in
 *                  PDB files
 *        Version:  1.0
 *        Created:  02/06/2008 01:31:11 PM CET
 *       Revision:  none
 *       Compiler:  g++
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */

#include <cstring>
#include <cstdlib>
#include "array.h"
#include "myfunc.h"
#include "mypro.h"

#define V_NUMRES 100 // variant of the number of residues between different records in pdb sequences
//#define MAX_NUM_FILE 90000
// ChangeLog
/*****************************************************************************
 * ChangeLog 2008-02-05
 *    AA3To1 using the advanced mode, that is, search the 3-letter residue name
 *    in PDB3Char_AllRes
 *
 * ChangeLog 2008-02-06
 *    add the functionality: for standard pdbfile names, set the chain ID in the same format as that in
 *    file pdb_seqres.txt, otherwise, use the old format
 *    the title and sequence format (AA_SEQ or DNA_SEQ) from PDB files is obtained
 * ChangeLog 2008-04-10
 *    when there are more than one line in TITLE item in the pdb file, only the
 *    annotation of the first line is achieved
 * ChangeLog 2009-06-10 
 *    Add the option --idtype, so the format of the id in the output file can
 *    be set.
 *    idtype == 0: 0 -- pdbid_chainID, e.g. 9wga_A
 *    idtype == 1: 0 -- PDBID+ChainID, e.g. 9WGAA
 ****************************************************************************/

//bool isStdPDBFileName = false; /*whether the supplied pdbfile name is standarized pdbfile name, that is 
//pdb$pdbid.ent or $pdbid.pdb, 2008-02-06, Nanjiang*/


/*****************************************************************************
 * Valgrind checked 2008-02-06, no leaks are possible
 ****************************************************************************/

bool isPrintDNASeq = false;

void PrintHelp()
{
    fprintf(stdout,"Usage: getseqresseq [options] pdbfile [chainid-list]\n");
    fprintf(stdout,"  get the protein sequence from SEQRES records pdb file\n");
    fprintf(stdout,"  if chainid is missing or 'all', all chains in the PDB file will be print out\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  -l file      : list file in the format\n");
    fprintf(stdout,"               : pdbfile chain-id-list, one record per line\n");
    fprintf(stdout,"  --dna        : print the DNA sequence as well\n");
    fprintf(stdout,"  -w N         : number of residues per line in output, default = 70\n");
    fprintf(stdout,"  -o file      : print the result to file, default = stdout\n");
    fprintf(stdout,"  --idtype 0|1 : set the idtype, default = 0\n");
    fprintf(stdout,"               : 0 -- pdbid_chainID, e.g. 9wga_A, 1 -- PDBID+ChainID, e.g. 9WGAA\n");
    //fprintf(stdout,"  -s|--std  : the pdbfile name is standard, i.e. pdb$pdbid.ent or $pdbid.pdb\n");
    //fprintf(stdout,"            : then the chain ID is named as lowercase($pdbid)_$chainID\n");
    fprintf(stdout,"  -h|--help : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-02-06, updated 2009-06-10,, Nanjiang Shu\n");
}
int GetSEQRESSeq(const char* pdbfile, Chain *chains, char *chainIDList, int &numChain)/*{{{*/
{
    int modeAA3To1  = 1; /*use the advanced mode to convert 3-char residue name to 1-char residue name, search in PDB3Char_AllRes, 2008-02-05, Nanjiang*/
    int i;
    FILE* fpPDBFile = NULL;
    fpPDBFile = fopen( pdbfile, "r");
    checkfilestream(fpPDBFile, pdbfile, "r");
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char recordid[SIZE_RECORD_ID+1] = "";
    int numResList[NUM_CHAINS] = { 0 }; // record the numRes for each chain as in the record of SEQRES 
    Chain *pChain = NULL;
    
    char proteinTitle[SIZE_LINE_PDB+1] = ""; /*title of the protein in the pdb file*/
    char *pTitle = NULL;

    while((linesize = fgetline(fpPDBFile, line ,maxline))!= EOF)
    {	
        if(sscanf(line,"%6s",recordid) != 1)
        {
            fprintf(stderr,"Bad pdb file!\n");
            fprintf(stderr,"line=\n%s\n",line);
            return -1;
        }
        if(pTitle == NULL && strncmp(recordid, "TITLE", 5) == 0) /*retrieve the protein annotation*/
            /*2008-04-10, add pTitle==NULL, so that only the first line
             * annotation is achieved*/
        {
            my_strcpy(proteinTitle, line+6, 64);
            pTitle = strtrim(proteinTitle);
        }
        else if (pTitle == NULL && strcmp(proteinTitle, "") == 0 && strncmp(recordid, "COMPND", 6) == 0)
        {
            my_strcpy(proteinTitle, line+6, 64);
            pTitle = strtrim(proteinTitle);
        }
        else if(strncmp(recordid,"SEQRES", SIZE_RECORD_ID)== 0 ) { break;}
    }

    char chainID;
    Array2D <char> resName_2darray(13,SIZE_RES_NAME+1);
    char **resName    = resName_2darray.array2D;               //initialize pointer as NULL

    //printf("SEQRES record found\n");
    // get the information about chains and residues of chains.
    int cntChain = -1 ; // count only the amino acid chain
    int serNum   = 0;
    int numRes   = 0;
    int resBegin = 0;
    int resEnd   = 0;
    char chainIDFormer = '\0';//initialize chainID as '\0'
    for(i = 0 ; i < NUM_CHAINS ; i++) InitChain(&chains[i]);
    while(1)/*{{{*/
    {
        Scanf_SEQRES_Para(line,recordid,serNum,chainID,numRes,resName[0]);				
        //			if(strcmp(id,"1A34")==0 && chainIDList[cntChain] == 'B')
        //				printf("%s\t%s\n",id,resName[0]);
        if(strncmp(recordid,"SEQRES", SIZE_RECORD_ID)== 0 // it should be in the SEQRES record
                //	&& strlen(resName[0]) == 3 // chains for amino acid
                && numRes >= 1 ) // the length of chain should be larger than 1
        {
            if(chainID != chainIDFormer) // the beginning of a new chain
            {
                if(cntChain > -1 )// before changing to the new chain
                {	
                    if(resEnd != numResList[cntChain])  
                        // check if the real numRes is equal to the numRes in the record
                        // if so, set the numRes for the chains as the residues actually
                        // exist in SEQRES record
                    {					
                        numResList[cntChain] = resEnd;
                        chains[cntChain].numRes = resEnd;					
                    }
                    chains[cntChain].aaSeq[chains[cntChain].numRes] = '\0';
                    //if(IsDNASeq(chains[cntChain].aaSeq, chains[cntChain].numRes))
                    //{
                    //    DeleteChain(&chains[cntChain]);
                    //    InitChain(&chains[cntChain]);
                    //    cntChain -- ;
                    //}
                }
                cntChain ++; // add a new chain

                chainIDList[cntChain] = chainID;
                numResList[cntChain] = numRes;

                chains[cntChain].chainID = chainID;
                chains[cntChain].numRes = numRes;
                chains[cntChain].aaSeq = new char[numRes+V_NUMRES+1];// in case of the real numRes is larger than the number in the record

                resBegin = 0 ; // the beginning index of aaSeq
                chainIDFormer = chainID;
            }
            resEnd = Scanf_SEQRES_Seq(line,resBegin,resName);

            for(i = resBegin; i < resEnd ; i ++)
                chains[cntChain].aaSeq[i] = AA3To1(resName[i-resBegin], modeAA3To1 );
            resBegin = resEnd;
        }
        else if(strncmp(recordid,"SEQRES", SIZE_RECORD_ID) != 0 )
        {
            if(cntChain > -1 )// before change to the new chain
            {							// check if the real numOfRes is 
                // equal to the numRes in the record
                if(resEnd < numResList[cntChain])  
                    // check if the real numOfRes is 
                {			   		// equal to the numRes in the record
                    numResList[cntChain] = resEnd;
                    chains[cntChain].numRes = resEnd;					
                }
                chains[cntChain].aaSeq[chains[cntChain].numRes] = '\0';
                //if(IsDNASeq(chains[cntChain].aaSeq,chains[cntChain].numRes))
                //{
                //    DeleteChain(&chains[cntChain]);
                //    InitChain(&chains[cntChain]);
                //    cntChain -- ;
                //}
            }
            cntChain ++;
            break;
        }
        if( (linesize = fgetline(fpPDBFile,line,maxline)) == EOF) break;
    }/*}}}*/

    fclose(fpPDBFile);
    chainIDList[cntChain] = '\0';
    numChain = cntChain;

    int sizeTitle = strlen(pTitle);
    for(i = 0 ; i < numChain; i ++)/*copy the title to each chain*/
    {
        pChain = &(chains[i]);
        pChain->title = new char[sizeTitle+1];
        my_strcpy(pChain->title, pTitle, sizeTitle);
        if(IsDNASeq(pChain->aaSeq, pChain->numRes))
        {
            pChain->seqtype = DNA_SEQ;
        }
        else 
        {
            pChain->seqtype = AA_SEQ;
        }
    }

    return numChain;
}
/*}}}*/
//int ReadInListFile(const char *listfile, char **pdbfilelist, char** chainidslist, int beginIndex)[>{{{<]
//{
//    int status = 0;
//    FILE* fp = NULL;
//    fp = fopen( listfile, "r");
//    int linesize;
//    int maxline = 300;
//    Array1D <char> line_1darray(maxline+1);
//    char *line = line_1darray.array1D;
//    int cntPDBFile  = beginIndex;
//    int linenumber = 0;
//    while((linesize = fgetline(fp, line ,maxline))!= EOF)
//    {
//        pdbfilelist[cntPDBFile] = new char[linesize+1];
//        chainidslist[cntPDBFile] = new char[linesize+1];
//        if(sscanf(line, "%s %s", pdbfilelist[cntPDBFile], chainidslist[cntPDBFile]) != 2)
//        {
//            fprintf(stderr, "listfile error at line %d\n", linenumber+1);
//            status = -1;
//            break;
//        }
//        linenumber ++;
//        cntPDBFile ++;
//    }
//    fclose(fp);

//    if(status != 0) return -1;
//    else return cntPDBFile;

//}
//[>}}}<]

int GetSelectedSEQRESSeq(const char *pdbfile, const char* chainids, FILE *fpout, int idtype = 0, int lineLength = 70)/*{{{*/
/*****************************************************************************
 * Get the amino acid sequences with specified chainids from SEQRES records
 * lineLength is the number of amino acid per line in output
 ****************************************************************************/
{
    int i = 0;
    Array1D <Chain> chains_1darray(NUM_CHAINS);
    Chain* chains = chains_1darray.array1D;

    for(i = 0 ; i < NUM_CHAINS; i++) InitChain(&chains[i]);
    char rtname[MAX_PATH+1] = "";
    char id[MAX_PATH+1] = "";
    char pdbid[SIZE_PDBID+1] = "";
    int numChain;
    char chainIDList[NUM_CHAINS+1] = "";
    bool isPrintAllChain = false;
    if(strcasecmp(chainids,"all") == 0 || strcasecmp(chainids,"") == 0)
        isPrintAllChain = true;

    numChain = GetSEQRESSeq(pdbfile, chains, chainIDList, numChain);
    for(i = 0; i < numChain ; i++)
    {
        if((IsInCharSet(chainIDList[i], chainids) || numChain <= 1 || isPrintAllChain) && (isPrintDNASeq || chains[i].seqtype == AA_SEQ))
        {                                                      
            rootname(pdbfile, rtname);
            if (StdPDBFileName2PDBID(rtname, pdbid) != NULL)
            {
                if (idtype == 0)
                {
                    my_strlwr(pdbid);
                    sprintf(id,"%s_%c", pdbid, chainIDList[i]);
                }
                else
                {
                    my_strupr(pdbid);
                    sprintf(id,"%s%c", pdbid, chainIDList[i]);
                }
            }
            else
            {
                my_strcpy(id, rtname, MAX_PATH);
            }
            fprintf(fpout, ">%s  mol:%s  length:%d  chainID:%c %s\n",id, 
                    (chains[i].seqtype == AA_SEQ)? "protein" : "nucleic", 
                    chains[i].numRes, chains[i].chainID, chains[i].title);
            WriteFastaSeq(chains[i].aaSeq, fpout, 0, 0x7FFFFFFF, lineLength);
        }
    }

    //free memory of chains
    for(i = 0 ; i < numChain; i++)
    {
        DeleteChain(&chains[i]);
    }
    return numChain;
}
/*}}}*/

int main(int argc, char** argv)/*{{{*/
{

    if(argc < 2)
    {
        PrintHelp();
        return -1;
    }

    int idtype = 0; /*idtype, 0 -- pdbid_chainID, e.g. 9wga_A, 1 -- PDBID+ChainID, e.g. 9WGAA, 2009-06-10*/
    char pdbfile[MAX_PATH+1] = "";
    char listfile[MAX_PATH+1] = "";
    char outfile[MAX_PATH+1] = "";
    char chainids[NUM_CHAINS+1] = "";
    int  lineLength = 70;
    int i = 1;
    const char *optionList[] = { "-h", "--help", "-l", "-o", "-w", "--line-length" , "--dna", "--idtype"};
    char numOption = sizeof(optionList)/sizeof(char*);

    bool isNonOptionArg = false;

    while (i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            {
                PrintHelp();
                return 0;
            }
            else if(strcmp(argv[i], "-l") == 0)
            {
                my_strcpy(listfile,argv[i+1], MAX_PATH);
                i +=2;
            }
            else if(strcmp(argv[i], "-o") == 0)
            {
                my_strcpy(outfile,argv[i+1], MAX_PATH);
                i +=2;
            }
            else if(strcmp(argv[i], "--idtype") == 0)
            {
                idtype = atoi(argv[i+1]);
                i +=2;
            }
            else if(strcmp(argv[i], "--dna") == 0)
            {
                isPrintDNASeq = true;
                i ++;
            }
            else if(strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--line-length") == 0)
            {
                lineLength = atoi(argv[i+1]);
                if(lineLength <= 0)
                {
                    fprintf(stderr, "Error! Line length can not be <= 0, reset to 70\n");
                    lineLength = 70;
                }
                i +=2;
            }
            else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            }
            else
            {
                fprintf(stderr,"Error! Wrong argument for option: '%s'\n", argv[i]);
                return -1;
            }
        }
        else //non-option argument
        {
            my_strcpy(pdbfile,argv[i], MAX_PATH);
            if(i+1 >= argc)
                break;
            else if(argv[i+1][0] == '-')
            {
                if(LinearSearch_String((const char*)argv[i+1], optionList, numOption) == -1)
                {
                    fprintf(stderr,"Error! Wrong argument for option: '%s'\n", argv[i+1]);
                    return -1;
                }
                else
                {
                    isNonOptionArg = false;
                    i ++;
                    continue;
                }
            }
            else
            {
                my_strcpy(chainids, argv[i+1], NUM_CHAINS);
                i += 2;
            }
        }
    }/*}}}*/

    FILE *fpout = NULL;
    if(strcmp(outfile,"") == 0)
    {
        fpout = stdout;
    }
    else
    {
        fpout = fopen(outfile,"w");
        checkfilestream(fpout, outfile, "w");
    }

    if(strcmp(pdbfile,"") != 0)
    {
        GetSelectedSEQRESSeq(pdbfile, chainids, fpout, idtype, lineLength);
    }
    else if(strcmp(listfile,"") != 0)
    {
        int linesize;
        int maxline = 300;
        Array1D <char> line_1darray(maxline+1);
        char *line = line_1darray.array1D;
        FILE *fpList;
        fpList = fopen(listfile,"r");
        checkfilestream(fpList,listfile,"r");

        while((linesize = fgetline(fpList, line, maxline)) != EOF)
        {
            if(linesize <= 0) continue;
            strcpy(pdbfile,"");
            strcpy(chainids,"");

            sscanf(line,"%s %s", pdbfile, chainids);
            GetSelectedSEQRESSeq(pdbfile, chainids, fpout, idtype, lineLength);
        }
        fclose(fpList);
    }
    else
    {
        fprintf(stderr,"Error! neither pdbfile and listfile is set!\n");
        return -1;
    }
    
    if(fpout != stdout && fpout != NULL) fclose(fpout);

    return EXIT_SUCCESS;
}
/*}}}*/

