/*
 * =====================================================================================
 *        Filename:  createseqmap.cpp
 *     Description:  create the indexing map between the sequence number and 
 *                   and the resSeq in PDB ATOM record
 *         Version:  1.0
 *         Created:  06/13/2006 03:44:25 PM CEST
 *        Revision:  none
 *        Compiler:  gcc
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */
/*****************************************************************************
 * valgrind checked, no leaks are possible, 2008-02-05, Nanjiang Shu
 ****************************************************************************/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
//#include <sys/stat.h>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"

int gapOpen  = -7;
int gapExt   = -1;
int match    = 5;
int misMatch = -10;
int MIN_MATCH_LENGTH = 3;
bool isPrintAlignment = false;
#define LOCAL_INIT_SHAPE '-'

#undef MAX_SIZE_ID 
#define MAX_SIZE_ID 100

#undef MAX_NUM_CHAIN
#define MAX_NUM_CHAIN 100

#undef V_NUMRES
#define V_NUMRES   200

//ChangeLog/*{{{*/
/*****************************************************************************
 * ChangeLog 2007-06-17
 * adding shape shatring, water accessibility, dssp secondary structure state
 * to seqmap file
 *
 * !! alignment left priority or right priority
 * seq-atm alignment, SS at 57 and 58, S57 is mutated at ATOM record 
 * but the alignment is left priority, so S57 in resseq-seq aligned to S58 in atom-seq
 * 
 * 661 : alignment between resseq seq and atom seq for id : 1H4XA
 * GapOpen      = -7
 * GapExtension = -1
 * identity     = 94.02%
 * similarity   = 94.02%
 * gapPercent   = 5.98%
 * identity of shorter seq= 100.00%
 * seq-A(117) <--> atm-A(110)
 * 
 *              123456789012345678901234567890123456789012345678901234567890   60
 * seq-A :    1 MAFQLEMVTRETVVIRLFGELDHHAVEQIRAKISTAIFQGAVTTIIWNFERLSFMDSSGV   60
 *               |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||
 * atm-A :    1 -AFQLEMVTRETVVIRLFGELDHHAVEQIRAKISTAIFQGAVTTIIWNFERLSFMDS-GV   58
 * 
 * 
 *              123456789012345678901234567890123456789012345678901234567  117
 * seq-A :   61 GLVLGRMRELEAVAGRTILLNPSPTMRKVFQFSGLGPWMMDATEEEAIDRVRGIVNG  117
 *              ||||||||||||||||||||||||||||||||||||||||||||||||||||     
 * atm-A :   59 GLVLGRMRELEAVAGRTILLNPSPTMRKVFQFSGLGPWMMDATEEEAIDRVR-----  110
 * 
 * ***************************************************
 * 661 : alignment between resseq seq and dssp seq for id : 1H4XA
 * GapOpen      = -7
 * GapExtension = -1
 * identity     = 94.02%
 * similarity   = 94.02%
 * gapPercent   = 5.13%
 * identity of shorter seq= 99.10%
 * seq-A(117) <--> dssp-A(111)
 * 
 *              123456789012345678901234567890123456789012345678901234567890   60
 * seq-A :    1 MAFQLEMVTRETVVIRLFGELDHHAVEQIRAKISTAIFQGAVTTIIWNFERLSFMDSSGV   60
 *               ||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||
 * dssp-A:    1 -AFQLEMVTRETVVIRLFGELDHHAVEQIRAKISTAIFQGAVTTIIWNFERLSFMDXSGV   59
 * 
 *              123456789012345678901234567890123456789012345678901234567  117
 * seq-A :   61 GLVLGRMRELEAVAGRTILLNPSPTMRKVFQFSGLGPWMMDATEEEAIDRVRGIVNG  117
 *              ||||||||||||||||||||||||||||||||||||||||||||||||||||     
 * dssp-A:   60 GLVLGRMRELEAVAGRTILLNPSPTMRKVFQFSGLGPWMMDATEEEAIDRVR-----  111
 *
 * some residues in ATOM record does not contain coordinations for eought
 * atoms, and DSSP will hence exclude such residues. for example, 1RYQA, 
 * first 5 residues, are all HIS. ATOM record keeps -5, -4 ,-3, -2, -1
 * dssp keeps -4,-3,-2. sequence alignment has problem to match them
 *
 * ***************************************************
 * ChangeLog 2007-10-11 
 * for scop domain sequence, id might be in two types
 * d3sdha_ : pdbid is 3sdh, chainid is A, single chain
 * d1a0i_1 : pdbid is 1a0i, chainid is null, part 1 of the chain
 * e1a0h.1a: pdbid is 1a0h, chainid is A, part one of the domain, which is chain a
 *
 * these three sequences were tested on 2007-10-19 in the folder $WORKDIR
 * ./createseqmap d3sdha_ --idtype 1 --aapath $CASIODATA3/wk/scop54/scop54aa/ -d ./
 * ./createseqmap d1a0i_1 --idtype 1 --aapath $CASIODATA3/wk/scop54/scop54aa/ -d ./
 * ./createseqmap e1a0h.1a --idtype 1 --aapath $CASIODATA3/wk/scop54/scop54aa/ -d ./
 * 
 * so pdbid = scopid[1:4]
 * if scopid[0] == 'd'
 * {
 *    chainid = scopid[5]
 * }
 * else if (scopid[0] == 'e')
 * {
 *    chainid = scopid[7]
 * }
 * 
 * ***************************************************
 * ChangeLog 2007-11-15
 *    Change the local INIT_SHAPE to '-'
 * 
 * ChangeLog 2008-02-05
 *    the functionality dealing with scop genetic domain sequences added.
 *    after SCOP 1.63, domain sequences with id start with 'g' contain
 *    multi-chains. To create seqmaps of these sequences, a new function
 *    CreateSeqMap_SCOPg was added. 
 *    each record should be one line, with
 *    ids chainIDs
 *    e.g. 
 *    g6r1x.1 ABCD
 *    meaning scop domain sequence g6r1x.1, this domain sequence contains four PDB chains
 *
 *    The algorithm is 
 *
 *    1. match all coor_chains (e.g. all A, B, C and D chains) to SCOP chains
 *       to retrieve waterAcc and secStruc information
 *
 *    2. align the target sequence to all coor_chains listed in the chainIDs
 *    (e.g all A, B, C and D chains), for aligned residue positions, set the
 *    waterAcc and secStruc as the matched positions of coor_chain. for non
 *    aligned residues, set to initialized values.
 *
 * ChangeLog 2008-07-06
 *    add the option --pdbaatype 0 | 1 
 *      0 -- 1bmf_A.aa
 *      1 -- 1BMFA.aa
 * ChangeLog 2008-07-07 
 *     assert(idtCnt / dssp_chain.numRes >= 0.9 )  changed to warning
 *     originally, it is used to check whether all residues with valid resSer
 *     in the PDB file should also exist in the DSSP file, so that the idtCnt
 *     should be exactly the same as dssp_chain.numRes. However, due to the
 *     different conversion method of the non-standard( e.g. mutated amino
 *     acid) residues, the difference between idtCnt and dssp_chain.numRes
 *     may be large, e.g. for the PDB chain 2IU4A
 * ChangeLog 2010-01-15
 *     1. bug fixed when creating seqmap files for SCOP domain sequences with the
 *     ID start with "g". In that case, this scop domain sequence will be
 *     aligned to all chains in the PDB file. However, some unwanted occasional
 *     matches will overwritten the previous correct matched. To solve this
 *     bug, the program will first check the quality of the alignment. With the
 *     hypothesis that the part of sequence used in the SCOP "g" starting
 *     domain sequences should be at least 10 amino acid long, any alignments
 *     without 10 gapless matches are occasional matches and should be ignored
 *     Also the value of  "scop_chain_cnt" should be checked to see if it is
 *     the same as "scop_chain.numRes"
 *
 *     2. The function for auto detect idtype is added
 * ChangeLog 2010-01-16 
 *     The options --idtype and --pdbaatype are removed. the idtype is
 *     auto detected. For the pdbaatype (or the naming scheme of the amino acid
 *     sequence files), there is only one format, that is the filename of the
 *     amino acid sequence is $id.aa
 * ChangeLog 2010-01-26
 *     For scop ids, if it starts with g, the script scopfasta2idlist.awk will
 *     generate the idlist from the annotation line, therefore, the idlist will
 *     be supplied.
 * ChangeLog 2011-03-15
 *     reset the default output for alignFile and logFile
 *     
 ****************************************************************************//*}}}*/

void PrintHelp()
{
    fprintf(stdout, "Usage: createseqmap [options] ids | -l idListFile\n");
    fprintf(stdout, "  Create the seqmap file given pdb chain id or scopid\n");
    fprintf(stdout, "  make sure that the PDB file, DSSP file and Shape String files are available\n");
    fprintf(stdout, "  idtype is auto detected\n");
    fprintf(stdout, "Format of the idListFile: one record per line\n");
    fprintf(stdout, "  g1ko6.1 AB\n");
    fprintf(stdout, "  d1jjda_\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "  -d | --outpath  <dir>  : output directory, default = ./\n");
    fprintf(stdout, "  -G | --gapOpen  <int>  : alignment parameter, gap open,      default = -7\n");
    fprintf(stdout, "  -E | --gapExt   <int>  : alignment parameter, gap extension, default = -1\n");
    fprintf(stdout, "  -M | --match    <int>  : alignment parameter, match score,   default =  5\n");
    fprintf(stdout, "  -S | --misMatch <int>  : alignment parameter, gap open,      default = 10\n");
    fprintf(stdout, "  --log           <file> : log file, default = /dev/stdout\n");
    fprintf(stdout, "  --align         <file> : output the alignment to file, default = /dev/null\n");
    fprintf(stdout, "  --ext           <file> : seqmap file extension, default = seqmap\n");
    //fprintf(stdout, "  --idtype   0|1|2     : 0 -- standardized id, 1 -- scop seq id\n");
    //fprintf(stdout, "                       : 2 -- scop seq id start with 'g' which includes multi-chains\n");
    //fprintf(stdout, "  --pdbaatype 0|1      : 0 -- e.g. 1bmf_A.aa , 1 -- e.g. 1BMFA.aa, default = 0\n");
    fprintf(stdout, "  --aapath     <dir>     : supply path for pdbaa files, default = ./\n");
    fprintf(stdout, "  --shapepath  <dir>     : supply path for shapestring files , default = ./\n");
    fprintf(stdout, "  -pa|--printalign       : print the alignment info to stdout\n");
    fprintf(stdout, "  -h |--help             : print this help message and exit\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Created 2006-06-13, updated 2011-03-15, Nanjiang Shu\n");

    fprintf(stdout, "Examples\n");
    fprintf(stdout, "   createseqmap d1a12a_ --aapath scopaa --shapepath pdbchainShapeString\n");
}

int GetSeq_SEQRES(const char* pdbfile, Chain *chains, char *chainIDList, int &numChain)/*{{{*/
/*****************************************************************************
 *  Get the aaSeq from the SEQRES record in the PDB file for chains specified
 *  in "chainIDList"
 *  when chainIDList == "", all chains existing in the SEQRES record are
 *  retrieved.
 *  the numChain return the actual number of chains retrieved
 *  Created 2010-01-15, updated 2010-01-15  
 ****************************************************************************/
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

    char tmpChainIDList[100] = "";
    Chain tmpChains[NUM_CHAINS];
    for (i= 0; i < NUM_CHAINS; i ++)
    {
        InitChain(&tmpChains[i]);
    }
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
    while(1)/*{{{*/
    {
        Scanf_SEQRES_Para(line,recordid,serNum,chainID,numRes,resName[0]);				
        //			if(strcmp(id,"1A34")==0 && tmpChainIDList[cntChain] == 'B')
        //				printf("%s\t%s\n",id,resName[0]);
        if(strncmp(recordid,"SEQRES", SIZE_RECORD_ID)== 0 // it should be in the SEQRES record
                //	&& strlen(resName[0]) == 3 // tmpChains for amino acid
                && numRes >= 1 ) // the length of chain should be larger than 1
        {
            if(chainID != chainIDFormer) // the beginning of a new chain
            {
                if(cntChain > -1 )// before changing to the new chain
                {	
                    if(resEnd != numResList[cntChain])  
                        // check if the real numRes is equal to the numRes in the record
                        // if so, set the numRes for the tmpChains as the residues actually
                        // exist in SEQRES record
                    {					
                        numResList[cntChain] = resEnd;
                        tmpChains[cntChain].numRes = resEnd;					
                    }
                    tmpChains[cntChain].aaSeq[tmpChains[cntChain].numRes] = '\0';
                    //if(IsDNASeq(tmpChains[cntChain].aaSeq, tmpChains[cntChain].numRes))
                    //{
                    //    DeleteChain(&tmpChains[cntChain]);
                    //    InitChain(&tmpChains[cntChain]);
                    //    cntChain -- ;
                    //}
                }
                cntChain ++; // add a new chain

                tmpChainIDList[cntChain] = chainID;
                numResList[cntChain] = numRes;

                tmpChains[cntChain].chainID = chainID;
                tmpChains[cntChain].numRes = numRes;
                if (tmpChains[cntChain].aaSeq == NULL)
                {
                    tmpChains[cntChain].aaSeq = new char[numRes+V_NUMRES+1];// in case of the real numRes is larger than the number in the record
                }

                resBegin = 0 ; // the beginning index of aaSeq
                chainIDFormer = chainID;
            }
            resEnd = Scanf_SEQRES_Seq(line,resBegin,resName);

            for(i = resBegin; i < resEnd ; i ++)
            {
                tmpChains[cntChain].aaSeq[i] = AA3To1(resName[i-resBegin], modeAA3To1 );
            }
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
                    tmpChains[cntChain].numRes = resEnd;					
                }
                tmpChains[cntChain].aaSeq[tmpChains[cntChain].numRes] = '\0';
                //if(IsDNASeq(tmpChains[cntChain].aaSeq,tmpChains[cntChain].numRes))
                //{
                //    DeleteChain(&tmpChains[cntChain]);
                //    InitChain(&tmpChains[cntChain]);
                //    cntChain -- ;
                //}
            }
            cntChain ++;
            break;
        }
        if( (linesize = fgetline(fpPDBFile,line,maxline)) == EOF) break;
    }/*}}}*/

    fclose(fpPDBFile);
    tmpChainIDList[cntChain] = '\0';
    int tmpNumChain = cntChain;

    int sizeTitle = strlen(pTitle);
    for(i = 0 ; i < tmpNumChain; i ++)/*copy the title to each chain*/
    {
        pChain = &(tmpChains[i]);
        if (pChain->title  == NULL)
        {
            pChain->title = new char[sizeTitle+1];
        }
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

    /*copy the aa sequences read from the PDB file to chains*/
    if (strcmp(chainIDList , "")== 0)
    {   /*copy all chains*/
        for (i = 0; i< tmpNumChain; i ++)
        {
            if (chains[i].title == NULL)
            {
                chains[i].title = new char [strlen(tmpChains[i].title)+1];
                strcpy(chains[i].title, "");
            }
            if (chains[i].aaSeq == NULL)
            {
                chains[i].aaSeq = new char [tmpChains[i].numRes+1];
                strcpy(chains[i].aaSeq, "");
            }
            CopyChain(&chains[i], &tmpChains[i]);
        }
        my_strcpy(chainIDList, tmpChainIDList, tmpNumChain);
        numChain = tmpNumChain;
    }
    else 
    {
        numChain = strlen(chainIDList);
        for (i = 0; i < numChain; i++)
        {
            int idx = Char2Digit(chainIDList[i], tmpChainIDList, tmpNumChain);
            if (idx != -1)
            {
                if (chains[i].title == NULL)
                {
                    chains[i].title = new char [strlen(tmpChains[idx].title)+1];
                    strcpy(chains[i].title, "");
                }
                if (chains[i].aaSeq == NULL)
                {
                    chains[i].aaSeq = new char [tmpChains[idx].numRes+1];
                    strcpy(chains[i].aaSeq , "" );
                }
                CopyChain(&chains[i], &tmpChains[idx]);
            }
            else 
            {
                fprintf(stderr,"Error! chainID '%c' does not exist in the PDB file '%s'\n", chainIDList[i], pdbfile);
            }
        }
    }
    
    /*clean allocated memory*/
    for (i = 0; i < NUM_CHAINS; i ++)
    {
        DeleteChain(&tmpChains[i]);
    }

    return numChain;
}
/*}}}*/
int GetPDBChainIDList(const char* pdbfile, char *chainIDList)/*{{{*/
{
    FILE* fpPDBFile = NULL;
    fpPDBFile = fopen( pdbfile, "r");
    checkfilestream(fpPDBFile, pdbfile, "r");
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char recordid[SIZE_RECORD_ID+1] = "";

    set <char> chainIDList_set;
    set <char> :: iterator css;
    
    
    bool isSEQRESArea = false;
    bool isEndSEQRESArea = false;

    while((linesize = fgetline(fpPDBFile, line ,maxline))!= EOF)
    {	
        if(sscanf(line,"%6s",recordid) != 1)
        {
            fprintf(stderr,"Bad pdb file!\n");
            fprintf(stderr,"line=\n%s\n",line);
            return -1;
        }
        if(strncmp(recordid,"SEQRES", SIZE_RECORD_ID)== 0 ) 
        {
            chainIDList_set.insert(line[11]);
            isSEQRESArea = true;

        }
        else if (isSEQRESArea)
        {
            isEndSEQRESArea = true;
        }
        if (isEndSEQRESArea)
        {
            break;
        }
    }
    int numChain = chainIDList_set.size();
    Set2Array(chainIDList_set.begin(), chainIDList_set.end(), chainIDList);
    chainIDList[numChain] = '\0';
    fclose(fpPDBFile);
    return numChain;
}
/*}}}*/
int DetectIDType(char *id)/*{{{*/
/*****************************************************************************
 *  detect the type of id
 *  idtype = 0: five-letter pdb chain id
 *  idtype = 1: scop domain id derived from single chain
 *  idtype = 2: scop domain id (starts with "g") derived from multiple chains
 *  idtype = 3: scop domain id (starts with "e")
 ****************************************************************************/
{
    int idtype =0;
    if (strlen (id) == 5)
    {
        idtype = 0;
    }
    else if (strlen (id) == 7)
    {
        if (id[0] == 'g')
        {
            idtype = 2;
        }
        else
        {
            idtype = 1;
        }
    }
    else if (strlen(id) == 8 && id [0] == 'e')
    {
        idtype = 3;
    }
    else 
    {
        idtype = -1; /*unidentified idtype*/
    }
    return idtype;
}/*}}}*/
int DetectChainID(char *id, char* chainIDList, int idtype )/*{{{*/
{
    int numChain = 1;
    if (idtype == 0)
    {
        chainIDList[0] = id[4];
    }
    else if (idtype == 1)
    { /*start with "d"*/
       chainIDList[0] = toupper(id[5]);
    }
    else if (idtype == 3)
    {  /*start with "e"*/
       chainIDList[0] = toupper(id[7]);
    }

    if (idtype == 2 || chainIDList[0] == '_')
    {
        char pdbid[SIZE_PDBID+1]= "";
        char pdbfilepath[MAX_PATH+1] = "";
        my_strcpy(pdbid, id+1, SIZE_PDBID);
        if(GetPDBFilePath(pdbid, pdbfilepath) != NULL)
        {
            numChain = GetPDBChainIDList(pdbfilepath, chainIDList);
        }
        else
        {
            fprintf(stderr,"Can not open pdbfile %s when processing the id %s\n",pdbfilepath, id);
            numChain  = -1;
        }
    }
    return numChain;
}/*}}}*/

int GetLongestGaplessMatch (int *alignRel, int alignLength) /*2010-01-15*//*{{{*/
/*****************************************************************************
 * Get the longest gapless match given an alignment
 ****************************************************************************/
{
    int i = 0;
    int max = 0;
    while (i < alignLength)
    {
        if (alignRel [i] == IDT)
        {
            int j = 0;
            while(alignRel[i+j] == IDT)
            {
                j ++;
            }
            if (j > max)
            {
                max = j;
            }
            i += j;
        }
        else
        {
            i += 1;
        }
    }

    return max;
}/*}}}*/
int CreateSeqMap(char* id, char chainID, const char* outpath, FILE * fpAlign, FILE *fpLog, int cnt ,  int idtype, const char *pdbaapath, const char *shapestringpath, char *seqmapfile)/*{{{*/
/*****************************************************************************
 * output to seqmapfile
 ****************************************************************************/
{

#ifdef DEBUG_ID
    printf("id=%s, idtype =%d\n", id, idtype);
#endif

    int j;

    char pdbaafilepath[MAX_PATH+1] = "";
    char pdbfilepath[MAX_PATH+1] = "";
    char shapestringfilepath[MAX_PATH+1] = "";
    char dsspfilepath[MAX_PATH+1] = "";


    Chain seqres_chain;
    Chain coor_chain;
    Chain dssp_chain;

    //Chain *pChain;
    InitChain(&seqres_chain);
    InitChain(&coor_chain);
    InitChain(&dssp_chain);
    seqres_chain.aaSeq    = new char [LONGEST_SEQ+1];
    seqres_chain.resSer   = new int [LONGEST_SEQ];
    seqres_chain.resICode = new char [LONGEST_SEQ+1];
    seqres_chain.waterAcc = new int [LONGEST_SEQ+1];
    seqres_chain.shString = new char [LONGEST_SEQ+1];
    seqres_chain.secStruc = new char [LONGEST_SEQ+1];

    coor_chain.aaSeq    = new char [LONGEST_SEQ+1];
    coor_chain.resSer   = new int [LONGEST_SEQ];
    coor_chain.resICode = new char [LONGEST_SEQ+1];

    dssp_chain.aaSeq    = new char [LONGEST_SEQ+1];
    dssp_chain.resSer   = new int [LONGEST_SEQ];
    dssp_chain.resICode = new char [LONGEST_SEQ+1];
    dssp_chain.waterAcc = new int [LONGEST_SEQ];
    dssp_chain.secStruc = new char [LONGEST_SEQ+1] ; 

    int lengthWholeChain = 0;
    Array1D <char> shString_1darray(LONGEST_SEQ+1);
    char *shString = shString_1darray.array1D;
    Array1D <char> aaSeqWholeChain_1darray(LONGEST_SEQ+1);
    char *aaSeqWholeChain = aaSeqWholeChain_1darray.array1D;

    int seq_type;
    char pdbid[SIZE_PDBID+1] = "";
    /*get pdbid*/
    if (idtype == 0)
    {
        my_strcpy(pdbid, id, SIZE_PDBID);
        my_strlwr(pdbid);
    }
    else 
    {
        my_strcpy(pdbid, id+1, SIZE_PDBID);
        my_strlwr(pdbid);
    }

    my_strcpy(seqres_chain.pdbid,pdbid,SIZE_PDBID);
    my_strcpy(coor_chain.pdbid,pdbid,SIZE_PDBID);
    my_strcpy(dssp_chain.pdbid,pdbid,SIZE_PDBID);
    seqres_chain.chainID = chainID;
    coor_chain.chainID   = chainID;
    dssp_chain.chainID   = chainID;

    bool isFileMissing = false; // check if the necessary files are prepared
    bool isShapeFileMissing = false;   // 2007-07-12, shape string missing is not that serious


    seq_type = UNKNOWN_SEQ_TYPE;

    sprintf(pdbaafilepath, "%s/%s.aa", pdbaapath, id);
    seqres_chain.numRes = ReadSeq_FASTA(pdbaafilepath, seqres_chain.aaSeq, &seq_type);
    my_strupr(seqres_chain.aaSeq);

    if (seqres_chain.numRes == READ_FILE_ERROR)
    {
        fprintf(stderr,"Can not open pdbaa file %s\n",pdbaafilepath);
        isFileMissing = true;
    }
    else if(seq_type == DNA_SEQ)
    {
        fprintf(stderr, "%d \t %s is a DNA sequence, neglect\n", cnt+1, id);
        isFileMissing = true;
    }

    if(GetPDBFilePath(pdbid,pdbfilepath) == NULL)/*{{{*/
    {
        fprintf(stderr,"Error! Can not open the pdbfile %s\n",pdbfilepath);
        isFileMissing = true;
    }/*}}}*/

    char stdid[SIZE_CHAIN_ID+1] = ""; /*standardized chain identifier, 5 characters*/
    char chainIDList[50+1] = "";
    sprintf(chainIDList,"%c", chainID);
    sprintf(stdid, "%s%c", pdbid, chainID);
    StdID(stdid);

    if(idtype != 0)/*{{{*/ /*read in the whole chain if the id is scopid*/
    {
        Chain tmpChain ;
        InitChain(&tmpChain);
        int tmpNumChain = 0;
        tmpNumChain = GetSeq_SEQRES(pdbfilepath, &tmpChain, chainIDList, tmpNumChain);
        if (tmpNumChain == -1)
        {
            fprintf(stderr,"Error! Can not read the chain '%c' from the PDB file %s\n",chainIDList[0], pdbfilepath);
            isFileMissing = true;
        }
        else
        {
            my_strcpy(aaSeqWholeChain, tmpChain.aaSeq, tmpChain.numRes);
        }
        lengthWholeChain = strlen(aaSeqWholeChain) ;
        DeleteChain(&tmpChain);
    }/*}}}*/

    GetSeq_ATOM(pdbfilepath,&coor_chain, chainIDList, false, true, LONGEST_SEQ, fpLog);
    /*GetSeq_ATOM return the number of chains read in*/
    if(IsDNASeq(coor_chain.aaSeq))
    {
        fprintf(stderr,"id %s is DNA sequence\n",id);
        isFileMissing = true;
    }

    if(GetShapeStringFilePath(stdid,shapestringfilepath, shapestringpath) != NULL)/*{{{*/
    {
        ReadSeq_FASTA(shapestringfilepath, shString, &seq_type);
        if(idtype == 0)
        {
            my_strcpy(seqres_chain.shString , shString, LONGEST_SEQ);
        }
        else if (idtype == 1)
        {
            // align aaSeqWholeChain and aaSeq of the scop sequence/*{{{*/
            //char *pch = strstr(aaSeqWholeChain, seqres_chain.aaSeq);
            //int shift = pch - aaSeqWholeChain;
            //if (pch != NULL)
            //{
            //    my_strcpy(seqres_chain.shString, shString+shift, seqres_chain.numRes);
            //}
            //else 
            //{
            //    fprintf(stderr, "Error! %s: scop sequence does not found in pdbaa sequence\n", id);
            //    return -1;
            //}
            char title1[20] = "seq-";
            char title2[20] = "wholech";
            int length1 = seqres_chain.numRes;
            int length2 = strlen(aaSeqWholeChain);
            Array1D <char>  alignXstr_1darray(length1 + length2 + 2); 
            Array1D <char>  alignYstr_1darray(length1 + length2 + 2);
            Array1D <int>   alignRel_1darray (length1 + length2 + 1);
            char *alignXstr = alignXstr_1darray.array1D;
            char *alignYstr = alignYstr_1darray.array1D;
            int  *alignRel  = alignRel_1darray.array1D;
            int alignLength = 0;
            if (fpAlign){
                fprintf(fpAlign,"***************************************************\n\n");
                fprintf(fpAlign,"%d : alignment between whole chain seq and scop seq for id : %s\n",cnt,id);
            }
            alignLength = SeqMatchAlign_Protein( seqres_chain.aaSeq, aaSeqWholeChain,title1,title2, alignXstr,alignYstr,alignRel,fpAlign, isPrintAlignment);
            if (fpAlign){
                fflush(fpAlign);
            }
            int seq_chain_cnt= 0 ;
            int whole_seq_cnt = 0 ;
            int cntIDT = 0;
            for(j = 0 ; j < alignLength ; j ++)
            {
                if(alignRel[j] == IDT)
                {
                    seqres_chain.shString[seq_chain_cnt]   = shString[whole_seq_cnt];
                    seq_chain_cnt ++;
                    whole_seq_cnt ++;
                    cntIDT ++;
                }
                else if(alignRel[j] == GAP)
                {
                    if(alignYstr[j] == CHAR_INDEL)
                    {
                        seqres_chain.shString[seq_chain_cnt] = LOCAL_INIT_SHAPE ; //set shString to '-', if this residue does not exist in whole_seq
                        seq_chain_cnt ++;
                    }
                    else if(alignXstr[j] == CHAR_INDEL) // if there are more residues in ATOM record than in SEQRES record, e.g. in 1EJG
                    {
                        whole_seq_cnt ++;
                        //	fprintf(fpLog,"false align, the SEQRES missing residue\n");
                        //	exit(0);
                    }
                }
                else if(alignRel[j] == MIS)
                {
                    seqres_chain.shString[seq_chain_cnt] = shString[whole_seq_cnt];
                    seq_chain_cnt ++;
                    whole_seq_cnt ++;
                }
            }/*}}}*/
            seqres_chain.shString[seqres_chain.numRes] = '\0';
            if (double (cntIDT) / min(length1 , length2) < 0.9)
            {
                fprintf(stderr, "Warning! %s : scop seq matched to seqres seq with less than < 0.9 identity\n", id);
            }
        }
    }
    else
    {
        fprintf(stderr,"Error! Can not open the shapestringfile %s\n",shapestringfilepath);
        isShapeFileMissing = true;
    }/*}}}*/
    if(GetDSSPFilePath(pdbid,dsspfilepath) != NULL)/*{{{*/
    {
        dssp_chain.numRes = GetDSSPChain(stdid, &dssp_chain, dsspfilepath);
        for(j = 0 ; j < dssp_chain.numRes; j ++)
        {
            if(islower(dssp_chain.aaSeq[j]))
                dssp_chain.aaSeq[j] = 'C';
        }
    }
    else
    {
        fprintf(stderr,"Can not open dsspfile %s\n",dsspfilepath);
        isFileMissing = true;
    }/*}}}*/

    if(!isFileMissing && ! isShapeFileMissing) /*shape string can not be missing either, 2008-02-05, Nanjiang*/
    {
        //if(isShapeFileMissing)
        //{
        //    for(j = 0 ; j < seqres_chain.numRes; j ++)
        //    {
        //        seqres_chain.shString[j] = LOCAL_INIT_SHAPE;
        //    }
        //    seqres_chain.shString[seqres_chain.numRes] = '\0';
        //}
        
        //alignment
        int  alignLength = 0;
        char title1[20] = "";
        char title2[20] = "";
        char title3[20] = "";
        int length1 = seqres_chain.numRes;
        int length2 = coor_chain.numRes;
        int length3 = dssp_chain.numRes;
        Array1D <char>  alignXstr_1darray(length1 + max(length2, length3) + 2); 
        Array1D <char>  alignYstr_1darray(length1 + max(length2, length3) + 2);
        Array1D <int>   alignRel_1darray (length1 + max(length2, length3) + 1);
        char *alignXstr = alignXstr_1darray.array1D;
        char *alignYstr = alignYstr_1darray.array1D;
        int  *alignRel  = alignRel_1darray.array1D;

        //Array1D <int>  resSer_from_dsspchain_1darray(seqres_chain.numRes);// used for checking consistency of the alignment between seq-atm and seq-dssp
        //int *resSer_from_dsspchain = resSer_from_dsspchain_1darray.array1D;

        sprintf(title1, "seq-%c", seqres_chain.chainID);
        sprintf(title2, "atm-%c", coor_chain.chainID);
        sprintf(title3, "dssp-%c", dssp_chain.chainID);

        // align resseq_seq and atom_seq/*{{{*/
        if(fpAlign){
            fprintf(fpAlign,"***************************************************\n\n");
            fprintf(fpAlign,"%d : alignment between resseq seq and atom seq for id : %s\n",cnt,id);
        }
        //fprintf(stdout,"***************************************************\n\n");
        //fprintf(stdout,"%d : alignment between resseq seq and atom seq for id : %s\n",cnt,id);
        alignLength = SeqMatchAlign_Protein( seqres_chain.aaSeq, coor_chain.aaSeq,title1,title2, alignXstr,alignYstr,alignRel,fpAlign, isPrintAlignment);
        if(fpAlign){
            fflush(fpAlign);
        }
        int seq_chain_cnt = 0 ;
        int coor_chain_cnt = 0 ;
        for(j = 0 ; j < alignLength ; j ++)
        {
            if(alignRel[j] == IDT)
            {
                seqres_chain.resSer[seq_chain_cnt]   = coor_chain.resSer[coor_chain_cnt];
                seqres_chain.resICode[seq_chain_cnt] = coor_chain.resICode[coor_chain_cnt];
                seq_chain_cnt ++;
                coor_chain_cnt ++;
            }
            else if(alignRel[j] == GAP)
            {
                if(alignYstr[j] == CHAR_INDEL)
                {
                    seqres_chain.resSer[seq_chain_cnt] = INIT_RESSEQ ; //set resSer to INIT_RESSEQ, if this residue does not exist in ATOM record
                    seqres_chain.resICode[seq_chain_cnt] = ' ';
                    seq_chain_cnt ++;
                }
                else if(alignXstr[j] == CHAR_INDEL) // if there are more residues in ATOM record than in SEQRES record, e.g. in 1EJG
                {
                    coor_chain_cnt ++;
                    //	fprintf(fpLog,"false align, the SEQRES missing residue\n");
                    //	exit(0);
                }
            }
            else if(alignRel[j] == MIS)
            {
                seqres_chain.resSer[seq_chain_cnt]   = coor_chain.resSer[coor_chain_cnt];
                seqres_chain.resICode[seq_chain_cnt] = coor_chain.resICode[coor_chain_cnt];
                seq_chain_cnt ++;
                coor_chain_cnt ++;
            }

        }/*}}}*/

        ////align resseq_seq and dssp_seq[>{{{<]
        //fprintf(fpAlign,"***************************************************\n\n");
        //fprintf(fpAlign,"%d : alignment between resseq seq and dssp seq for id : %s\n",cnt,id);
        //fprintf(stdout,"***************************************************\n\n");
        //fprintf(stdout,"%d : alignment between resseq seq and dssp seq for id : %s\n",cnt,id);
        //alignLength = SeqMatchAlign_Protein( seqres_chain.aaSeq, dssp_chain.aaSeq,title1,title3, alignXstr,alignYstr,alignRel,fpAlign);
        //fflush(fpAlign);
        //seq_chain_cnt = 0;
        //int dssp_chain_cnt = 0;
        //for(j = 0 ; j < alignLength ; j ++)
        //{
        //    if(alignRel[j] == IDT)
        //    {
        //        seqres_chain.waterAcc[seq_chain_cnt] = dssp_chain.waterAcc[dssp_chain_cnt];
        //        seqres_chain.secStruc[seq_chain_cnt] = dssp_chain.secStruc[dssp_chain_cnt];

        //        resSer_from_dsspchain[seq_chain_cnt] = dssp_chain.resSer[dssp_chain_cnt];
        //        seq_chain_cnt ++;
        //        dssp_chain_cnt ++;
        //    }
        //    else if(alignRel[j] == GAP)
        //    {
        //        if(alignYstr[j] == CHAR_INDEL)
        //        {
        //            seqres_chain.waterAcc[seq_chain_cnt] = INIT_WATERACC;
        //            seqres_chain.secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;

        //            resSer_from_dsspchain[seq_chain_cnt] = INIT_RESSEQ;
        //            seq_chain_cnt ++;
        //        }
        //        else if(alignXstr[j] == CHAR_INDEL) // if there are more residues in ATOM record than in SEQRES record, e.g. in 1EJG
        //        {
        //            dssp_chain_cnt ++;
        //        //  fprintf(fpLog,"false align, the SEQRES missing residue\n");
        //        }
        //    }
        //    else if(alignRel[j] == MIS)
        //    {
        //        seqres_chain.waterAcc[seq_chain_cnt] = dssp_chain.waterAcc[dssp_chain_cnt];
        //        seqres_chain.secStruc[seq_chain_cnt] = dssp_chain.secStruc[dssp_chain_cnt];

        //        resSer_from_dsspchain[seq_chain_cnt] = dssp_chain.resSer[dssp_chain_cnt];
        //        seq_chain_cnt ++;
        //        dssp_chain_cnt ++;
        //    }
        //}[>}}}<]

        int dssp_chain_cnt = 0;/*{{{*/
        seq_chain_cnt  = 0;
        int idtcnt  = 0; // count for matched residue 
        while(1)// match dssp_seqres.resSer to seqres_chain.resSer, which is derived from coor_chain.resSer by sequence alignment. assume all resSer in DSSP exist in coor.resSer
        {
            if(seq_chain_cnt >= seqres_chain.numRes || dssp_chain_cnt >= dssp_chain.numRes) break; //bug fixed, && changed to ||, 2007-06-20

            if(seqres_chain.resSer[seq_chain_cnt] == INIT_RESSEQ)
            {
                seqres_chain.waterAcc[seq_chain_cnt] = INIT_WATERACC;
                seqres_chain.secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;
                seq_chain_cnt ++;
            }
            else
            {
                if(seqres_chain.resSer[seq_chain_cnt] == dssp_chain.resSer[dssp_chain_cnt] && seqres_chain.resICode[seq_chain_cnt] == dssp_chain.resICode[dssp_chain_cnt])
                {
                    seqres_chain.waterAcc[seq_chain_cnt] = dssp_chain.waterAcc[dssp_chain_cnt]; 
                    seqres_chain.secStruc[seq_chain_cnt] = dssp_chain.secStruc[dssp_chain_cnt];
                    dssp_chain_cnt ++;
                    seq_chain_cnt ++;
                    idtcnt ++;
                }
                else if(seqres_chain.resSer[seq_chain_cnt] > dssp_chain.resSer[dssp_chain_cnt])
                {
                    dssp_chain_cnt ++;
                }
                else
                {
                    seqres_chain.waterAcc[seq_chain_cnt] = INIT_WATERACC;
                    seqres_chain.secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;
                    seq_chain_cnt ++;
                }
            }
        }/*}}}*/
        while(seq_chain_cnt < seqres_chain.numRes)
        {
            seqres_chain.waterAcc[seq_chain_cnt] = INIT_WATERACC;
            seqres_chain.secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;
            seq_chain_cnt ++;
        }

#ifdef DEBUG_READDSSP
            fprintf(stdout,"#DSSP Chain: %s%c\n", dssp_chain.pdbid,dssp_chain.chainID);
            for(j = 0; j < dssp_chain.numRes; j ++)
            {
                fprintf(stdout,"%c %3d %c %c\n", dssp_chain.aaSeq[j], dssp_chain.resSer[j], dssp_chain.resICode[j], dssp_chain.chainID);
            }
            fprintf(stdout,"\n");

            fprintf(stdout,"#seqres Chain: %s%c\n", seqres_chain.pdbid, seqres_chain.chainID);
            for(j = 0; j < seqres_chain.numRes; j ++)
            {
                fprintf(stdout,"%c %3d %c %c\n", seqres_chain.aaSeq[j], seqres_chain.resSer[j], seqres_chain.resICode[j], seqres_chain.chainID);
            }
            fprintf(stdout,"\n");
#endif
        if (idtype == 0)
        {
            if(idtcnt / double(dssp_chain.numRes) < 0.9)
            {
                fprintf(stderr, "Warning! for chain %s%c dssp_chain vs seqres_chain,  idtCnt  = %d, idtCnt/dssp_chain.numRes <0.9\n", dssp_chain.pdbid, dssp_chain.chainID, idtcnt);
                fprintf(stderr, "dssp_chain_cnt = %d, dssp_chain.numRes   = %d\n", dssp_chain_cnt, dssp_chain.numRes);
                fprintf(stderr, "seq_chain_cnt  = %d, seqres_chain.numRes = %d\n", seq_chain_cnt, seqres_chain.numRes);
                //fflush(stdout);
                //assert( idtcnt / double(dssp_chain.numRes) >= 0.9) ;
            }
        }
        
        ////for(j = 0; j < seqres_chain.numRes; j++)// check the consistency of the alignment of[>{{{<]
        //    // seqres_chain.aaSeq and coor_chain.aaSeq
        //    // seqres_chain.aaSeq and dssp_chain.aaSeq 
        ////{
        //    //if(seqres_chain.resSer[j] != resSer_from_dsspchain[j] &&
        //            //seqres_chain.resSer[j] != INIT_RESSEQ &&
        //            //resSer_from_dsspchain[j] != INIT_RESSEQ)
        //    //{
        //        //fprintf(stdout, "seqresResSer[%d]= %d, resSerFromDSSP[%d]= %d",j, seqres_chain.resSer[j],j, resSer_from_dsspchain[j] );
        //        //assert(seqres_chain.resSer[j] == resSer_from_dsspchain[j]) ; // debug
        //    //}
        //    //NOTE, the sequence in ATOM record of dssp and might not be the
        //    //same, for example, in 1A5T, the last PRO330 in ATOM record does
        //    //not exist in the dssp file
        //    //
        //    //it could be that the beginning residue in DSSP does not exist in
        //    //ATOM reocord in PDB file, for example, 1A62, MSE1 is annotated as
        //    //HETATM, while in DSSP, it is listed as the first residue of the
        //    //sequence
        //    //
        //    //when there are two consecutive same residue in the sequence and
        //    //one of them is mutated, for example, in 1H4XA, SER57 is mutated
        //    //to be SEP in ATOM record, in the aligment, 
        //    //
        ////}[>}}}<]

        // output the seqmap file/*{{{*/
        int length = seqres_chain.numRes;
        FILE *fpout = fopen(seqmapfile,"w");
        checkfilestream(fpout, seqmapfile, "w");
        int k;
        fprintf(fpout,"# ID: %s\n",id);
        fprintf(fpout,"# sequence length: %d\n",length);
        fprintf(fpout,"# AA       -- 1 letter amino acid\n");
        fprintf(fpout,"# seqidx   -- sequence index for the sequence in pdb_seqres.txt, starting from 1\n");
        fprintf(fpout,"# resser   -- index of the residue in the pdb file, \n");
        fprintf(fpout,"#             %d means the residue exists in the sequence, but not exists in ATOM record\n", INIT_RESSEQ);
        fprintf(fpout,"# icode    -- pdb insertion code for the residue, '%c' means icode is blank\n", NULL_ICODE);
        fprintf(fpout,"# shape    -- 8-state shape symbol, '%c' means shape symbol is not availble\n", UNKNOWN_SHAPE);
        fprintf(fpout,"# waterAcc -- water accessibility from DSSP, '%d' means water accessibility is not availble\n", INIT_WATERACC);
        fprintf(fpout,"# dsspSec  -- secondary structure symbol in DSSP, '%c' means secondary structure is not availble\n", UNKNOWN_DSSP_SEC);
        fprintf(fpout,"#             '%c' means random coil, which is blank in DSSP\n", DSSP_SEC_RANDOM);
        fprintf(fpout,"%2s %6s %6s %5s %5s %7s %7s\n","AA", "seqidx", "resser","icode", "shape", "waterAcc", "dsspSec");
        for( k = 0 ; k < length ; k ++)
        {
            char iCode ;
            char secStruc ;
            if((iCode = seqres_chain.resICode[k]) == ' ')
                iCode = NULL_ICODE;
            if((secStruc = seqres_chain.secStruc[k]) == ' ')
                secStruc = DSSP_SEC_RANDOM;
            fprintf(fpout,"%-2c %6d %6d %5c %5c %7d %7c\n",seqres_chain.aaSeq[k], k+1, seqres_chain.resSer[k], iCode, seqres_chain.shString[k], seqres_chain.waterAcc[k], secStruc);
        }
        fclose(fpout);/*}}}*/

    }
    //free memory
    DeleteChain(&seqres_chain);
    DeleteChain(&coor_chain);
    DeleteChain(&dssp_chain);


    return 0;
}/*}}}*/
int CreateSeqMap_SCOPg(char* scopid, char* chainIDList, const char* outpath, FILE * fpAlign, FILE *fpLog, int cnt ,  int idtype, const char *pdbaapath, const char *shapestringpath, char *seqmapfile)/*{{{*/
/*****************************************************************************
 * create seqmap file for SCOP genetic domain sequences, with ids starting with
 * 'g'. These domain sequences compose multi-chains
 * 2008-02-05, Nanjiang Shu
 * output to seqmapfile
 ****************************************************************************/
{

#ifdef DEBUG_ID
    printf("scopid=%s, idtype =%d\n", scopid, idtype);
#endif

    int i,j;

    char scopseqfilepath[MAX_PATH+1] = "";
    char pdbfilepath[MAX_PATH+1] = "";
    char shapestringfilepath[MAX_PATH+1] = "";
    char dsspfilepath[MAX_PATH+1] = "";

    int numChain = strlen(chainIDList); /* number of chains in this scop genetic domain. There are multi-chains in the scop genetic domain,2008-02-05, Nanjiang*/

    /*get the PDBID*/
    char pdbid[SIZE_PDBID+1] = "";
    my_strcpy(pdbid, scopid+1, SIZE_PDBID);
    my_strlwr(pdbid);

    Chain scop_chain;  /*scop_chain got from Astral Databse http://astral.berkeley.edu */
    Array1D <Chain> seqres_chain_1darray(numChain);
    Array1D <Chain> coor_chain_1darray(numChain);
    Array1D <Chain> dssp_chain_1darray(numChain);
    Chain *seqres_chain = seqres_chain_1darray.array1D;  /*pdb chains generated from SEQRES records in PDB files*/
    Chain *coor_chain = coor_chain_1darray.array1D;        /*pdb chains generated from ATOM records in PDB files*/
    Chain *dssp_chain = dssp_chain_1darray.array1D;        /*dssp chains generated from DSSP files*/

    Chain *pChain;
    InitChain (&scop_chain);
    for(i = 0 ; i < numChain; i ++)
    {
        InitChain(&seqres_chain[i]);
        InitChain(&coor_chain[i]);
        InitChain(&dssp_chain[i]);
    }

    /*allocate memory for chains*/
    pChain = &scop_chain; 
    pChain->aaSeq    = new char [LONGEST_SEQ+1];
    pChain->resSer   = new int [LONGEST_SEQ];
    pChain->resICode = new char [LONGEST_SEQ+1];
    pChain->waterAcc = new int [LONGEST_SEQ+1];
    pChain->shString = new char [LONGEST_SEQ+1];
    pChain->secStruc = new char [LONGEST_SEQ+1];
    for(i = 0 ; i < numChain; i ++)/*{{{*/
    {
        pChain = &(seqres_chain[i]);
        pChain->aaSeq    = new char [LONGEST_SEQ+1];
        pChain->resSer   = new int [LONGEST_SEQ];
        pChain->resICode = new char [LONGEST_SEQ+1];
        pChain->shString = new char[LONGEST_SEQ+1];
        pChain->waterAcc = new int [LONGEST_SEQ];
        pChain->secStruc = new char [LONGEST_SEQ+1] ; 

        pChain = &(coor_chain[i]);
        pChain->aaSeq    = new char [LONGEST_SEQ+1];
        pChain->resSer   = new int [LONGEST_SEQ];
        pChain->resICode = new char [LONGEST_SEQ+1];

        pChain = &(dssp_chain[i]);
        pChain->aaSeq    = new char [LONGEST_SEQ+1];
        pChain->resSer   = new int [LONGEST_SEQ];
        pChain->resICode = new char [LONGEST_SEQ+1];
        pChain->waterAcc = new int [LONGEST_SEQ];
        pChain->secStruc = new char [LONGEST_SEQ+1] ; 
    }/*}}}*/

    int seq_type;
    char chainID = ' ';
    char stdid[SIZE_CHAIN_ID+1] = ""; /*standardized chain identifier, 5 characters*/

    for(i = 0; i < numChain; i ++)
    {
        chainID = chainIDList[i];
        my_strcpy(seqres_chain[i].pdbid,pdbid,SIZE_PDBID);
        my_strcpy(coor_chain[i].pdbid,pdbid,SIZE_PDBID);
        my_strcpy(dssp_chain[i].pdbid,pdbid,SIZE_PDBID);
        seqres_chain[i].chainID = chainID;
        coor_chain[i].chainID   = chainID;
        dssp_chain[i].chainID   = chainID;
    }

    bool isFileMissing = false; /*check if the necessary files are prepared*/
    bool isShapeFileMissing = false;   /* 2007-07-12, shape string missing is not that serious*/
    /*first get the SCOP domain sequences*//*{{{*/
    sprintf(scopseqfilepath, "%s/%s.aa", pdbaapath, scopid); /*idtype = 2*/
    scop_chain.numRes = ReadSeq_FASTA(scopseqfilepath, scop_chain.aaSeq, &seq_type);
    my_strupr(scop_chain.aaSeq);
    if (scop_chain.numRes == READ_FILE_ERROR)
    {
        fprintf(stderr,"Error! Can not open scop seq file %s\n",scopseqfilepath);
        isFileMissing = true;
    }
    else if(seq_type == DNA_SEQ)
    {
        fprintf(stderr, "%d \t %s is a DNA sequence, neglect\n", cnt+1, scopid);
        isFileMissing = true;
    }/*}}}*/

    int len = 0; /*temporary storing sequence length*/

    if(GetPDBFilePath(pdbid,pdbfilepath) == NULL)
    {
        fprintf(stderr,"Error! Can not open the PDB file %s\n",pdbfilepath);
        isFileMissing = true;
    }

    /*read in aaSeq for all chains*/
    int tmpNumChain = 0;
    tmpNumChain = GetSeq_SEQRES(pdbfilepath, seqres_chain, chainIDList, tmpNumChain);
    if (tmpNumChain != numChain)
    {
        fprintf(stderr,"Error! The number of chains read from the pdbfile %s = %d is not equal to the number of chains (%d) need to be read in\n", pdbfilepath, tmpNumChain, numChain);
    }/*2010-01-15*/

    for(i = 0 ; i < numChain; i ++) /*{{{*/ /*get seqres_chain, coor_chain and dssp_chain*/
    {
        sprintf(stdid, "%s%c", pdbid, chainIDList[i]);
        StdID(stdid);
        pChain = &seqres_chain[i] ;  /*get the seqres_chains for all chains included by the SCOP chain*/

        if(GetShapeStringFilePath(stdid,shapestringfilepath, shapestringpath) != NULL)/*{{{*/
        {
            len = ReadSeq_FASTA(shapestringfilepath, pChain->shString, &seq_type); /*here pChain = &seqres_chain[i]*/
            if (len != pChain->numRes)
            {
                fprintf(stderr,"Error! %s seqres_seq (%d) and shapestring (%d) not the same length\n", stdid, pChain->numRes, len);
                assert (len == pChain->numRes);
            }
        }
        else
        {
            fprintf(stderr,"Error! Can not open the shape string file %s\n",shapestringfilepath);
            isShapeFileMissing = true;
        }/*}}}*/

        pChain = &(dssp_chain[i]); /*{{{*/
        if(GetDSSPFilePath(pdbid,dsspfilepath) != NULL)
        {
            pChain->numRes = GetDSSPChain(stdid, pChain, dsspfilepath);
            for(j = 0 ; j < pChain->numRes; j ++)
            {
                if(islower(pChain->aaSeq[j]))
                    pChain->aaSeq[j] = 'C';
            }
        }
        else
        {
            fprintf(stderr,"Error! Can not open the dssp file %s\n",dsspfilepath);
            isFileMissing = true;
        }/*}}}*/

    }/*}}}*/

    GetSeq_ATOM(pdbfilepath,coor_chain, chainIDList, false, true, LONGEST_SEQ, fpLog);
    /*GetSeq_ATOM return the number of chains read in*/
    for(i = 0 ; i < numChain; i ++)
    {
        pChain = &(coor_chain[i]);
        if(IsDNASeq(pChain->aaSeq))
        {
            fprintf(stderr,"Error! id %s%c is a DNA sequence\n",pChain->pdbid, pChain->chainID);
            isFileMissing = true;
        }
    }

    if(!isFileMissing && !isShapeFileMissing)
    {
#ifdef DEBUG_READSEQ
        fprintf(stdout, "chainIDList = %s\n", chainIDList);
        for(i = 0 ; i < numChain; i ++)
        {
            fprintf(stdout, "seq-%c: %s\n", seqres_chain[i].chainID, seqres_chain[i].aaSeq);
            fprintf(stdout, "atm-%c: %s\n", coor_chain[i].chainID, coor_chain[i].aaSeq);
            fprintf(stdout,"\n");

        }
#endif
        //if(isShapeFileMissing)/*{{{*//*obsolete code*/
        //{
        //    for(j = 0 ; j < seqres_chain.numRes; j ++)
        //    {
        //        seqres_chain.shString[j] = LOCAL_INIT_SHAPE;
        //    }
        //    seqres_chain.shString[seqres_chain.numRes] = '\0';
        //}/*}}}*/


        /*first align coor_chain and dssp_chain to seqres_chain for each chain in chainIDList*/
        for(i = 0; i < numChain; i ++)/*{{{*/
        {
            //alignment
            int  alignLength = 0;
            char title1[20] = "";
            char title2[20] = "";
            char title3[20] = "";
            int length1 = seqres_chain[i].numRes;
            int length2 = coor_chain[i].numRes;
            int length3 = dssp_chain[i].numRes;
            Array1D <char>  alignXstr_1darray(length1 + max(length2, length3) + 2); 
            Array1D <char>  alignYstr_1darray(length1 + max(length2, length3) + 2);
            Array1D <int>   alignRel_1darray (length1 + max(length2, length3) + 1);
            char *alignXstr = alignXstr_1darray.array1D;
            char *alignYstr = alignYstr_1darray.array1D;
            int  *alignRel  = alignRel_1darray.array1D;

            //Array1D <int>  resSer_from_dsspchain_1darray(seqres_chain[i].numRes);// used for checking consistency of the alignment between seq-atm and seq-dssp
            //int *resSer_from_dsspchain = resSer_from_dsspchain_1darray.array1D;

            sprintf(title1, "seq-%c", seqres_chain[i].chainID);
            sprintf(title2, "atm-%c", coor_chain[i].chainID);
            sprintf(title3, "dssp-%c", dssp_chain[i].chainID);

            // align resseq_seq and atom_seq/*{{{*/
            if (fpAlign){
                fprintf(fpAlign,"***************************************************\n\n");
                fprintf(fpAlign,"%d : alignment between resseq seq and atom seq for scopid : %s\n",cnt,scopid);
            }
            //fprintf(stdout,"***************************************************\n\n");
            //fprintf(stdout,"%d : alignment between resseq seq and atom seq for scopid : %s, chainID = %c\n",cnt,scopid, seqres_chain[i].chainID);
            alignLength = SeqMatchAlign_Protein( seqres_chain[i].aaSeq, coor_chain[i].aaSeq,title1,title2, alignXstr,alignYstr,alignRel,fpAlign, isPrintAlignment);
            if(fpAlign){
                fflush(fpAlign);
            }
            int seq_chain_cnt = 0 ;
            int coor_chain_cnt = 0 ;
            for(j = 0 ; j < alignLength ; j ++)
            {
                if(alignRel[j] == IDT)
                {
                    seqres_chain[i].resSer[seq_chain_cnt]   = coor_chain[i].resSer[coor_chain_cnt];
                    seqres_chain[i].resICode[seq_chain_cnt] = coor_chain[i].resICode[coor_chain_cnt];
                    seq_chain_cnt ++;
                    coor_chain_cnt ++;
                }
                else if(alignRel[j] == GAP)
                {
                    if(alignYstr[j] == CHAR_INDEL)
                    {
                        seqres_chain[i].resSer[seq_chain_cnt] = INIT_RESSEQ ; //set resSer to INIT_RESSEQ, if this residue does not exist in ATOM record
                        seqres_chain[i].resICode[seq_chain_cnt] = ' ';
                        seq_chain_cnt ++;
                    }
                    else if(alignXstr[j] == CHAR_INDEL) // if there are more residues in ATOM record than in SEQRES record, e.g. in 1EJG
                    {
                        coor_chain_cnt ++;
                        //	fprintf(fpLog,"false align, the SEQRES missing residue\n");
                        //	exit(0);
                    }
                }
                else if(alignRel[j] == MIS)
                {
                    seqres_chain[i].resSer[seq_chain_cnt]   = coor_chain[i].resSer[coor_chain_cnt];
                    seqres_chain[i].resICode[seq_chain_cnt] = coor_chain[i].resICode[coor_chain_cnt];
                    seq_chain_cnt ++;
                    coor_chain_cnt ++;
                }
            }/*}}}*/

            int dssp_chain_cnt = 0;/*{{{*/
            seq_chain_cnt  = 0;
            int idtcnt  = 0; // count for matched residue 
            while(1)// match dssp_seqres.resSer to seqres_chain.resSer, which is derived from coor_chain.resSer by sequence alignment. assume all resSer in DSSP exit in coor.resSer
            {
                if(seq_chain_cnt >= seqres_chain[i].numRes || dssp_chain_cnt >= dssp_chain[i].numRes) break; //bug fixed, && changed to ||, 2007-06-20

                if(seqres_chain[i].resSer[seq_chain_cnt] == INIT_RESSEQ)
                {
                    seqres_chain[i].waterAcc[seq_chain_cnt] = INIT_WATERACC;
                    seqres_chain[i].secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;
                    seq_chain_cnt ++;
                }
                else
                {
                    if(seqres_chain[i].resSer[seq_chain_cnt] == dssp_chain[i].resSer[dssp_chain_cnt] && seqres_chain[i].resICode[seq_chain_cnt] == dssp_chain[i].resICode[dssp_chain_cnt])
                    {
                        seqres_chain[i].waterAcc[seq_chain_cnt] = dssp_chain[i].waterAcc[dssp_chain_cnt]; 
                        seqres_chain[i].secStruc[seq_chain_cnt] = dssp_chain[i].secStruc[dssp_chain_cnt];
                        dssp_chain_cnt ++;
                        seq_chain_cnt ++;
                        idtcnt ++;
                    }
                    else if(seqres_chain[i].resSer[seq_chain_cnt] > dssp_chain[i].resSer[dssp_chain_cnt])
                    {
                        dssp_chain_cnt ++;
                    }
                    else
                    {
                        seqres_chain[i].waterAcc[seq_chain_cnt] = INIT_WATERACC;
                        seqres_chain[i].secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;
                        seq_chain_cnt ++;
                    }
                }
            }/*}}}*/
            while(seq_chain_cnt < seqres_chain[i].numRes)
            {
                seqres_chain[i].waterAcc[seq_chain_cnt] = INIT_WATERACC;
                seqres_chain[i].secStruc[seq_chain_cnt] = UNKNOWN_DSSP_SEC;
                seq_chain_cnt ++;
            }

#ifdef DEBUG_READDSSP
            fprintf(stdout,"#DSSP Chain: %s%c\n", dssp_chain[i].pdbid,dssp_chain[i].chainID);
            for(j = 0; j < dssp_chain[i].numRes; j ++)
            {
                fprintf(stdout,"%c %3d %c %c\n", dssp_chain[i].aaSeq[j], dssp_chain[i].resSer[j], dssp_chain[i].resICode[j], dssp_chain[i].chainID);
            }
            fprintf(stdout,"\n");

            fprintf(stdout,"#seqres Chain: %s%c\n", seqres_chain[i].pdbid, seqres_chain[i].chainID);
            for(j = 0; j < seqres_chain[i].numRes; j ++)
            {
                fprintf(stdout,"%c %3d %c %c\n", seqres_chain[i].aaSeq[j], seqres_chain[i].resSer[j], seqres_chain[i].resICode[j], seqres_chain[i].chainID);
            }
            fprintf(stdout,"\n");
#endif

            if (idtype == 0)
            {
                if(idtcnt / double(dssp_chain[i].numRes) < 0.9)
                {
                    fprintf(stderr, "Warning! for chain %s%c dssp_chain vs seqres_chain,  idtCnt  = %d, idtCnt/dssp_chain.numRes <0.9\n", dssp_chain[i].pdbid, dssp_chain[i].chainID, idtcnt);
                    fprintf(stderr, "dssp_chain_cnt = %d, dssp_chain.numRes   = %d\n", dssp_chain_cnt, dssp_chain[i].numRes);
                    fprintf(stderr, "seq_chain_cnt  = %d, seqres_chain.numRes = %d\n", seq_chain_cnt, seqres_chain[i].numRes);
                    //fflush(stdout);
                    //assert( idtcnt / double(dssp_chain[i].numRes) >= 0.9) ;
                }
            }

            ////for(j = 0; j < seqres_chain.numRes; j++)// check the consistency of the alignment of[>{{{<]
            //    // seqres_chain.aaSeq and coor_chain.aaSeq
            //    // seqres_chain.aaSeq and dssp_chain.aaSeq 
            ////{
            //    //if(seqres_chain.resSer[j] != resSer_from_dsspchain[j] &&
            //            //seqres_chain.resSer[j] != INIT_RESSEQ &&
            //            //resSer_from_dsspchain[j] != INIT_RESSEQ)
            //    //{
            //        //fprintf(stdout, "seqresResSer[%d]= %d, resSerFromDSSP[%d]= %d",j, seqres_chain.resSer[j],j, resSer_from_dsspchain[j] );
            //        //assert(seqres_chain.resSer[j] == resSer_from_dsspchain[j]) ; // debug
            //    //}
            //    //NOTE, the sequence in ATOM record of dssp and might not be the
            //    //same, for example, in 1A5T, the last PRO330 in ATOM record does
            //    //not exist in the dssp file
            //    //
            //    //it could be that the beginning residue in DSSP does not exist in
            //    //ATOM reocord in PDB file, for example, 1A62, MSE1 is annotated as
            //    //HETATM, while in DSSP, it is listed as the first residue of the
            //    //sequence
            //    //
            //    //when there are two consecutive same residue in the sequence and
            //    //one of them is mutated, for example, in 1H4XA, SER57 is mutated
            //    //to be SEP in ATOM record, in the aligment, 
            //    //
            ////}[>}}}<]

        }/*}}}*/

        /*initialize the shString, waterAcc and secStruc in scop_chain*/
        pChain = &scop_chain;/*{{{*/
        for(j = 0 ; j < scop_chain.numRes; j ++)
        {
            pChain->resSer[j]   = INIT_RESSEQ ;
            pChain->resICode[j] = NULL_ICODE;
            pChain->waterAcc[j] = INIT_WATERACC;
            pChain->shString[j] = UNKNOWN_SHAPE;
            pChain->secStruc[j] = UNKNOWN_DSSP_SEC;
        }/*}}}*/

/*========================================================*/
        /*then align several seqres_chain to scop_chain*/
        int cntMatchedRes = 0; /*2010-01-15, count the number of matched residues for the scop domain sequence*/
        for(i = 0; i < numChain; i ++)/*{{{*/
        {
            //alignment
            int  alignLength = 0;
            char title1[20] = "";
            char title2[20] = "";
            int length1 = seqres_chain[i].numRes;
            int length2 = scop_chain.numRes;
            Array1D <char>  alignXstr_1darray(length1 + length2+ 2); 
            Array1D <char>  alignYstr_1darray(length1 + length2+ 2);  
            Array1D <int>   alignRel_1darray (length1 + length2+ 2);
            char *alignXstr = alignXstr_1darray.array1D;
            char *alignYstr = alignYstr_1darray.array1D;
            int  *alignRel  = alignRel_1darray.array1D;

            sprintf(title1, "resseq-%c", seqres_chain[i].chainID);
            sprintf(title2, "scop");

            // align resseq_seq[i] and scop_chain
            if (fpAlign){
                fprintf(fpAlign,"***************************************************\n\n");
                fprintf(fpAlign,"%d : alignment between resseq seq and scop seq for scopid : %s, chainID = %c\n",cnt,scopid, seqres_chain[i].chainID);
            }
            //fprintf(stdout,"***************************************************\n\n");
            //fprintf(stdout,"%d : alignment between resseq seq and scop seq for scopid : %s, chainID = %c\n",cnt,scopid, seqres_chain[i].chainID);
            alignLength = SeqMatchAlign_Protein( seqres_chain[i].aaSeq, scop_chain.aaSeq,title1,title2, alignXstr,alignYstr,alignRel,fpAlign, isPrintAlignment);
            if(fpAlign){
                fflush(fpAlign);
            }

            int longestGaplessMatch = 100;
            longestGaplessMatch = GetLongestGaplessMatch (alignRel, alignLength); /*2010-01-15*/

#ifdef DEBUG
            printf("iChain=%d, alignLength=%d, numChain= %d, longestGaplessMatch=%d\n",i, alignLength, numChain,longestGaplessMatch); //debug
#endif 

            if (longestGaplessMatch < MIN_MATCH_LENGTH)
            {
                continue;
            }
            int seq_chain_cnt = 0 ;
            int scop_chain_cnt = 0 ;
            int continuousMatchRes = 0; /*2010-01-16 , only those alignment with more than MIN_MATCH_LENGTH continuous residues are considered as real match*/
            j = 0;
            while(j  < alignLength)    /*there is no SIM for seqmatchalign_protein*/
            {
                if(alignRel[j] == IDT || alignRel[j] == MIS)
                {
                    int k = 0;
                    while((j+k < alignLength)&&(alignRel[j+k] != GAP))
                    {
                        k ++;
                    }
                    continuousMatchRes = k;

                    //printf("iChain=%d, alignLength=%d, numChain= %d, continuousMatchRes=%d\n",i, alignLength, numChain,continuousMatchRes); //debug

                    if (continuousMatchRes >= MIN_MATCH_LENGTH)
                    {
                        for (k = j; k < j + continuousMatchRes; k ++)
                        {
                            scop_chain.resSer[scop_chain_cnt]  = seqres_chain[i].resSer[seq_chain_cnt];
                            scop_chain.resICode[scop_chain_cnt]  = seqres_chain[i].resICode[seq_chain_cnt];
                            scop_chain.shString[scop_chain_cnt]  = seqres_chain[i].shString[seq_chain_cnt];
                            scop_chain.waterAcc[scop_chain_cnt]  = seqres_chain[i].waterAcc[seq_chain_cnt];
                            scop_chain.secStruc[scop_chain_cnt]  = seqres_chain[i].secStruc[seq_chain_cnt];
#ifdef DEBUG
                            fprintf(stdout,"scopAA vs seqresAA [%d] = %c %c    %c %c\n", k, scop_chain.aaSeq[scop_chain_cnt], seqres_chain[i].aaSeq[seq_chain_cnt], alignXstr[k], alignYstr[k]);
#endif
                            scop_chain_cnt ++;
                            seq_chain_cnt ++;
                            cntMatchedRes ++;
                        }
                    }
                    else
                    {
                        scop_chain_cnt += continuousMatchRes;
                        seq_chain_cnt += continuousMatchRes;
                    }
                    j += continuousMatchRes;
                }
                else // if(alignRel[j] == GAP)
                {
                    if(alignYstr[j] == CHAR_INDEL)
                    {
                        seq_chain_cnt ++;
                    }
                    else if(alignXstr[j] == CHAR_INDEL) // if there are more residues in ATOM record than in SEQRES record, e.g. in 1EJG
                    {
                        scop_chain_cnt ++;
                    }
                    j ++;
                }
            }
        }/*}}}*/

        if (cntMatchedRes != scop_chain.numRes)
        {
            fprintf(stdout,"Warning! cntMatchedRes (%d) != scop_chain.numRes (%d) \n", cntMatchedRes, scop_chain.numRes);
        }

        /* output the seqmap file *//*{{{*/
        int length = scop_chain.numRes;
        FILE *fpout = fopen(seqmapfile,"w");
        if (checkfilestream(fpout, seqmapfile, "w") != -1)
        {
            int k;
            fprintf(fpout,"# scopid: %s\n",scopid);
            fprintf(fpout,"# sequence length: %d\n",length);
            fprintf(fpout,"# AA       -- 1 letter amino acid\n");
            fprintf(fpout,"# seqidx   -- sequence index for the sequence in pdb_seqres.txt, starting from 1\n");
            fprintf(fpout,"# resser   -- index of the residue in the pdb file, \n");
            fprintf(fpout,"#             %d means the residue exists in the sequence, but not exists in ATOM record\n", INIT_RESSEQ);
            fprintf(fpout,"# icode    -- pdb insertion code for the residue, '%c' means icode is blank\n", NULL_ICODE);
            fprintf(fpout,"# shape    -- 8-state shape symbol, '%c' means shape symbol is not availble\n", UNKNOWN_SHAPE);
            fprintf(fpout,"# waterAcc -- water accessibility from DSSP, '%d' means water accessibility is not availble\n", INIT_WATERACC);
            fprintf(fpout,"# dsspSec  -- secondary structure symbol in DSSP, '%c' means secondary structure is not availble\n", UNKNOWN_DSSP_SEC);
            fprintf(fpout,"#             '%c' means random coil, which is blank in DSSP\n", DSSP_SEC_RANDOM);
            fprintf(fpout,"%2s %6s %6s %5s %5s %7s %7s\n","AA", "seqidx", "resser","icode", "shape", "waterAcc", "dsspSec");
            for( k = 0 ; k < length ; k ++)
            {
                char iCode ;
                char secStruc ;
                if((iCode = scop_chain.resICode[k]) == ' ')
                    iCode = NULL_ICODE;
                if((secStruc = scop_chain.secStruc[k]) == ' ')
                    secStruc = DSSP_SEC_RANDOM;
                fprintf(fpout,"%-2c %6d %6d %5c %5c %7d %7c\n",scop_chain.aaSeq[k], k+1, scop_chain.resSer[k], iCode, scop_chain.shString[k], scop_chain.waterAcc[k], secStruc);
            }
            fclose(fpout);
        }/*}}}*/
    }
    //free memory
    DeleteChain(&scop_chain);
    for(i = 0; i< numChain; i ++)
    {
        DeleteChain(&seqres_chain[i]);
        DeleteChain(&coor_chain[i]);
        DeleteChain(&dssp_chain[i]);
    }

    return 0;
}/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    if( argc < 2 )
    {
       fprintf(stdout, "too few arguments\n");
       PrintHelp();
       return -1;
    }

    int  i;
    char idListFile[MAX_PATH+1] = "";
    char ext[MAX_PATH+1]        = "seqmap";
    char outpath[MAX_PATH+1]    = "./";
    char alignFile[MAX_PATH+1]  = "null";
    char logFile[MAX_PATH+1]    = "/dev/stdout";
    char datadir[MAX_PATH+1]    = "";
    char id[MAX_SIZE_ID+1]    = "";    /*allocate more memory for id, since scopid can be > 8, 2007-10-11 */
    int idtype = 0; /*default idtype = 0, meaning standardized 5 character chain identifier*/
    //int pdbaatype  = 0; [>default = 0, the naming type of pdbaafile, 0 -- 1bmf_A.aa 1 -- 1BMFA.aa 2008-07-06  <]
    char pdbaapath[MAX_PATH+1] = "./";
    char shapestringpath[MAX_PATH+1] = "./";

    GetDataDir(datadir);
    //     sprintf(outpath,"%s/%s",datadir,"seqmap"); //default outpath = $DATADIR/seqmap

    vector <string> v_id;
    vector <string> v_chainIDList;
    vector <string> :: iterator ivs;

    bool isNonOptionArg = false;

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
            {
                PrintHelp();
                return 0;
            }
            else if(strcmp(argv[i],"-d") == 0 || strcmp(argv[i], "--outpath") == 0)
            {
                my_strcpy(outpath,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--align") == 0 )
            {
                my_strcpy(alignFile,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"-pa") == 0 || strcmp(argv[i],"--printalign") == 0 )
            {
                isPrintAlignment = true;
                i += 1;
            }
            else if(strcmp(argv[i],"--ext") == 0 )
            {
                my_strcpy(ext,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--gapOpen") == 0|| strcmp(argv[i], "-G") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, gapOpen, false)) == -1)
                    return -1;
            }
            else if(strcmp(argv[i],"--gapExt") == 0 || strcmp(argv[i], "-E" ) == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, gapExt, false)) == -1)
                    return -1;
            }
            else if(strcmp(argv[i],"--match") == 0|| strcmp(argv[i], "-M" ) == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, match, false)) == -1)
                    return -1;
            }
            else if(strcmp(argv[i],"--misMatch") == 0|| strcmp(argv[i], "-S" ) == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, misMatch, false)) == -1)
                    return -1;
            }
            //else if(strcmp(argv[i],"--idtype") == 0)
            //{
            //if( ( i = option_parser_numeric(argc, argv, i, idtype, true, 0, 3)) == -1)
            //return -1;
            //}
            //else if(strcmp(argv[i],"--pdbaatype") == 0)
            //{
                //if( ( i = option_parser_numeric(argc, argv, i, pdbaatype, true, 0, 1)) == -1)
                    //return -1;
            //}
            else if(strcmp(argv[i],"--log") == 0|| strcmp(argv[i], "-log" ) == 0)
            {
                my_strcpy(logFile,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--aapath") == 0)
            {
                my_strcpy(pdbaapath,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--shapepath") == 0)
            {
                my_strcpy(shapestringpath,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"-l") == 0)
            {
                my_strcpy(idListFile,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            }
            else
            {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        }
        else
        { 
            idtype = DetectIDType(argv[i]);
            if (idtype == 2) 
            {
                if (i+1 < argc && argv[i+1][0] != '-' )
                {
                    v_id.push_back(argv[i]);
                    v_chainIDList.push_back(argv[i+1]);
                    i += 2;
                }
                else
                {
                    v_id.push_back(argv[i]);
                    v_chainIDList.push_back("");
                    i ++;
                }
            }
            else if (idtype >= 0)
            {
                v_id.push_back(argv[i]);
                v_chainIDList.push_back("");
                i ++;
            }
            else 
            {
                fprintf(stderr,"Unrecognized ID = %s, please check! Exit ...", argv[i]);
                return EXIT_FAILURE;
            }
        }
    }/*}}}*/

    if( strcmp(idListFile,"") == 0 && v_id.size() <= 0)
    {
        fprintf(stderr, "Error! Neither idListFile nor ids set in the argument list\n");
        return -1;
    }

    VerifyFolder(outpath);

    FILE *fpLog = NULL;
    FILE *fpAlign = NULL;

    if (strcasecmp (logFile, "/dev/stdout") == 0 || strcasecmp(logFile, "stdout")==0) {
        fpLog = stdout;
    } else if (strcasecmp (logFile, "/dev/stderr") == 0 || strcasecmp(logFile, "stderr")==0) {
        fpLog = stderr;
    } else if (strcasecmp(logFile, "") != 0 && strcasecmp(logFile,"null") != 0){
        fpLog = fopen(logFile,"w");
        checkfilestream(fpLog, logFile,"w");
    }
    if (strcasecmp (alignFile, "/dev/stdout") == 0 || strcasecmp(alignFile, "stdout")==0) {
        fpAlign = stdout;
    } else if (strcasecmp (alignFile, "/dev/stderr") == 0 || strcasecmp(alignFile, "stderr")==0) {
        fpAlign = stderr;
    } else if (strcasecmp(alignFile,"") != 0 && strcasecmp(alignFile,"null") != 0){
        fpAlign = fopen(alignFile,"w");
        checkfilestream(fpAlign, alignFile,"w");
    }

    int cntRecord = 0;

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    if(strcmp(idListFile,"") !=0) // if idListFile set, input the ids in the idListFile to ids_set
    {
        FILE *fpin;
        fpin = fopen(idListFile, "r");
        checkfilestream(fpin, idListFile, "r");
        while((linesize = fgetline(fpin, line, maxline)) != EOF)
        {
            if (linesize <=0)
            {
                continue;
            }
            char *pch; 
            pch = strtok (line, WHITE_SPACE);
            int cntField=0;
            while (pch != NULL)
            {
                cntField ++;
                if (cntField == 1)
                {
                    v_id.push_back(pch);
                }
                else if (cntField == 2)
                {
                    v_chainIDList.push_back(pch);
                }
                pch = strtok(NULL, WHITE_SPACE);
            }
            if (cntField <= 0)
            {
                continue;
            }
            else if (cntField < 2)
            {
                v_chainIDList.push_back("");
            }
        }
        fclose(fpin);
    }

    //for(i = 0; i < int (v_id.size()); i ++)
    //{
        //fprintf(stderr, "%s %s\n", v_id[i].c_str(), v_chainIDList[i].c_str());
    //}

    for(i = 0; i < int (v_id.size()); i ++)
    {
        my_strcpy(id, v_id[i].c_str(), MAX_SIZE_ID);
        idtype = DetectIDType(id);

        if (idtype < 0)
        {
            fprintf(stderr,"Error! Unrecognized id %s\n", id);
            continue;
        }

        int numChain = 0;
        char chainIDList[200] = "";
        char seqmapfile[MAX_PATH+1] = "";
        if (v_chainIDList[i].compare("") == 0)
        {
            numChain = DetectChainID(id, chainIDList, idtype);
        }
        else
        {
            my_strcpy(chainIDList, v_chainIDList[i].c_str(), 100);
            numChain=strlen(chainIDList);
        }

        sprintf(seqmapfile,"%s/%s.%s",outpath, id, ext); 

        if (idtype == 2)
        {
            CreateSeqMap_SCOPg(id,chainIDList, outpath, fpAlign, fpLog, cntRecord, idtype, pdbaapath, shapestringpath, seqmapfile);
        }
        else
        {
            CreateSeqMap(id,chainIDList[0], outpath, fpAlign, fpLog, cntRecord,  idtype, pdbaapath, shapestringpath, seqmapfile);
        }

        fprintf(stdout,"%d\t %s seqmap output to %s\n",cntRecord+1, id, seqmapfile);
        cntRecord ++;
    }

    if (fpLog != NULL && fpLog != stderr && fpLog != stdout){
        fclose(fpLog);
    }
    if (fpAlign != NULL && fpAlign != stderr && fpAlign != stdout){
        fclose(fpAlign);
    }

    return EXIT_SUCCESS;
}
/*}}}*/
