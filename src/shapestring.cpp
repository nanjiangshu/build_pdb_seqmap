#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "myfunc.h"
#include "mypro.h"
#include "mytemplate.h"
#include "array.h"


#ifndef _strupr
#define strupr my_strupr
#define _strupr strupr
#endif 

// This code is under development not finished yet

struct PointDef 
{
    unit8 x;
    unit8 y;
    unit8 intensity;
};
int MAXSeqLen = 20000;
int printLineLength = 60; /*maximum length of each line to print out the shape string*/
bool isForceRunDSSP = false;
bool isUsingPDBAtom = false;
bool isChainIDCaseSensitive = true; /*by default the chainID is case sensitive*/

/*****************************************************************************
 * valgrind checked, not error. 2007-11-07, Nanjiang
 ****************************************************************************/
//ChangeLog/*{{{*/

/*****************************************************************************
 * ChangeLog 2007-10-21
 * fixed the program in reading sequence in SEQRES record, the trailing blanks
 * should be matched with "   ", but ""
 *
 * fixed the reference sequence for outputing the result, using caseq1 instead
 * of caseq2
 *
 * ChangeLog 2007-11-05
 * add the option for supplying the external amino acid sequence, by default,
 * the sequence is from SEQRES record
 *
 * ChangeLog reading in the binary data in rawdata
 *
 * ChangeLog 2007-11-07
 * Using the aigment function SeqMatchAlign_Protein, bug fixed finally
 *
 * ChangeLog 2007-11-09
 * using the torsion angles from DSSP by default, otherwise using
 * --forcedssp means force running dssp program to generate the dssp file
 *
 * calculate the shapestring after alignment, copying the phi psi angles
 * after the alignment. From torsion angles to shapestring, we need to know the
 * amino acid type also
 *
 * ChangeLog 2007-12-21
 * when dssp file name is supplied, pdbid should be generated as well
 *
 * ChangeLog 2008-07-07 
 *    add the option --nocase, an option to control the case sensitivity of the
 *    chainID, by default the chainID is case sensitive
 * ChangeLog 2008-08-30
 *    the default data path changed to $DATADIR/shapestring/rawdata instead of
 *    /misc/casiodata2/shapestring/rawdata
 *    The datatype of shapestring_rawdata_path changed from const char* to
 *    char*
 * ChangeLog 2009-12-21
 *    C$ChainID changed to C$ChainIDList, multiple chainIDs can be supplied in
 *    one line. Also, if $ChainIDList is empty, shape strings are created for
 *    all chains in the pdb file.
 *    ==under development, this functionality has not been finished yet
 * ChangeLog 2010-04-14 
 *    The datatype for PointDef changed to unit8 from int
 * 
 ****************************************************************************//*}}}*/
const char *AAs[] = /*{{{*/
{
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "PRO",
    "PHE",
    "MET",
    "LYS",
    "ARG",
    "HIS",
    "GLY",
    "SER",
    "THR",
    "CYS",
    "TYR",
    "ASN",
    "GLU",
    "TRP",
    "ASP",
    "GLN",
    "SUM",
    "A18"
};/*}}}*/

int sizeAAalpabet = 22;
char AAalphabet[40] = "";

void PrintHelp()
{
    fprintf(stdout,"Usage: shapesting [options] pdbfile [C$ChainIDList title aaSeqFile dsspfile] \n") ;
    fprintf(stdout,"  Use torsion angles from dssp, if dssp file does not exist, envoke the dssp program\n");
    fprintf(stdout,"  Note that the chainID is case sensitive by default, or use the option --nocase to disable it\n");
    fprintf(stdout,"  when a single chainID is supplied, 'title' is printed in the header line of the output file\n");
    fprintf(stdout,"  when multiple chainIDs are supplied, 'title$chainID' is printed.\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"options: \n");
    fprintf(stdout,"  -l <listfile> : list file for pdbfiles\n");
    fprintf(stdout,"  -o <outfile>  : output to the file, default is to stdout\n");
    fprintf(stdout,"  -a            : write the amino acid sequence above shape line\n");
    fprintf(stdout,"  --forcedssp   : force using dssp programs to generate dssp file\n");
    fprintf(stdout,"  --pdb         : calculate torsion angles from the ATOM record directly\n");
    fprintf(stdout,"  --rtxt        : read the text format rawdata, default readin the binary format\n");
    fprintf(stdout,"  --data <path> : path for the raw data, default=$DATADIR/shapestring/rawdata\n");
    fprintf(stdout,"  --max-seq int : setting the longest length of the sequence, default = 10000\n");
    fprintf(stdout,"  --linesize int: setting the linesize to print out shape string, default = 60\n");
    fprintf(stdout,"  --align  file : output the alignment to file, default not output the alignment\n");
    fprintf(stdout,"  --nocase      : set the chainID as case insensitive, default is case sensitive\n");
    fprintf(stdout,"  -h|--help     : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"filelist-name should be in the format: \n");
    fprintf(stdout,"filename CchainID [title]\n");
    fprintf(stdout,"*** one record per line\n");
    fprintf(stdout,"when title is not supplied, using pdbfilename as title\n");
    fprintf(stdout,"if aaSeqFile is not supplied, using the sequence from RESSEQ record\n");
    fprintf(stdout,"for exmaple\n");
    fprintf(stdout,"/pdbpath/pdb1abc.ent CA 1ABCA #Ca means using only chain a\n");
    fprintf(stdout,"/pdbpath/pdb1abd.ent CBCD 1ABDB #Cbcd means creating shape strings for chain b, c and d\n");
    fprintf(stdout,"/pdbpath/pdb1abd.ent C  1ABDB #C means creating shape strings for all chains\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Note:\n");
    fprintf(stdout,"  extract pdbid from pdbfile, pdbid=substr(rootname(pdbfile):3:6)\n");
    fprintf(stdout,"  search dssp file in $DATADIR/dssp\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Example:\n");
    fprintf(stdout,"    shapestring /data/1bmf.pdb CA 1BMF\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created 2006-01-01, updated 2010-04-14\n");
}

int GetSeq_SEQRES(const char* pdbfile, char chainID, char *seq)/*{{{*/
{
    FILE *fpPDBFile = NULL;
    fpPDBFile = fopen(pdbfile,"r");
    if(checkfilestream(fpPDBFile, pdbfile, "r") == -1)
    {
        return -1;
    }
    char recordid [SIZE_RECORD_ID+1]  ="";
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    int resBegin = 0;
    int resEnd = 0;
    int i;
    Array2D <char> resName_2darray(13,SIZE_RES_NAME+1);
    char **resName    = resName_2darray.array2D;               //initialize pointer as NULL

    while((linesize = fgetline(fpPDBFile, line, maxline)) != EOF)
    {
        sscanf(line,"%6s",recordid);
        if (  (strcmp(recordid,"SEQRES")==0) && (line[11]==chainID)  )  /*read in sequence from SEQRES record, 2007-10-21, nanjiang*/
        {
            resEnd = Scanf_SEQRES_Seq(line,resBegin,resName);

            for(i = resBegin; i < resEnd ; i ++)
                seq[i] = AA3To1(resName[i-resBegin]);
            resBegin = resEnd;

        }
        else if(strcmp(recordid, "ATOM") == 0)
        {
            break;
        }
    }
    fclose(fpPDBFile);
    seq[resEnd] = '\0';
    int seqLen = resEnd;
    return seqLen;
}
/*}}}*/

int GetParameter(const char *shapestring_rawdata_path, int8*** Shapestring, int *Align_Table, bool isReadTextRawdata = false)/*{{{*/
/*****************************************************************************
 * read in Align_Table
 * and initialize Shapestring[][][]
 ****************************************************************************/
{
    int ik, ij;
    int Ntemp;
    int nbeg;
    int nend;
    char openname[MAX_PATH+1] = "";
    FILE *fp;
    if (isReadTextRawdata)
    {
        for (ik=0; ik<20; ik++)
        {
            sprintf(openname,"%s/Table%s.txt",shapestring_rawdata_path, AAs[ik]);
            fp = fopen(openname,"r");
            if ( fp == NULL )
            {
                return 10;
            }
            while (fscanf(fp,"%d%d%d\n",&nbeg, &nend,&Ntemp ) != EOF)
            {
                if ( Ntemp >= 10 )
                {
                    Ntemp = Ntemp - 10;
                }
                Shapestring[ik][nbeg][nend] = int8 (Ntemp) ;
            }
            fclose(fp);
        }
    }
    else /*read in binary file*/ 
    {
        for (ik=0; ik<20; ik++)
        {
            sprintf(openname,"%s/Table%s.bin",shapestring_rawdata_path, AAs[ik]);
            fp = fopen(openname,"rb");
            if ( fp == NULL )
            {
                return 10;
            }
            int num = 0;
            size_t nread = 0; /*number of counts read in by function read()*/
            nread = fread(&num, sizeof(int), 1 , fp);
            if (nread != 1)
            {
                fprintf(stderr,"Read bindary file error! File:%s\n", openname);
            }
            if(num <0 || num > 5000000)
            {
                fprintf(stderr, "rawdata file %s corrupted\n", openname);
                assert(num >= 0 && num < 5000000);
            }
            int j =0 ;

            Array1D <PointDef> points_1darray(num);
            PointDef *points = points_1darray.array1D;
            nread = fread(points, sizeof(PointDef), num, fp);
            if (nread != size_t(num))
            {
                fprintf(stderr,"Read bindary file error! File:%s\n", openname);
            }
            fclose(fp);

            for(j = 0 ; j < num; j ++)
            {
                if ( points[j].intensity >= 10 )
                {
                    points[j].intensity -= 10;
                }
                Shapestring[ik][points[j].x][points[j].y] = int8 (points[j].intensity) ;
            }
        }
    }

    for (ik=0; ik<20; ik++)
    {
        for (ij=0; ij<20; ij++)
        {
            Align_Table[ik*20+ij] = -10;
            if (  ik == ij  )
            {
                Align_Table[ik*20+ij] =4;

            }
        }
    }
    Align_Table[400] = -7;//gap open
    Align_Table[401] = -2;//gap extension vertical, first sequence
    Align_Table[402] = -2;//gap extension horizontal, second sequence
    Align_Table[403] = -2;//unknown AAS

    return 0;
}/*}}}*/
double Calculate_torsion_angles_phi_psi(double *point1,double *point2,double *point3,double *point4)/*{{{*/
{
    double x10,y10,z10,x20,y20,z20,x30,y30,z30,cx1,cy1,cz1;
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,cx2,cy2,cz2;

    double costheta;
    double nx;
    double ny;
    double nz;
    double nx1;
    double ny1;
    double nz1;
    double xcrossd;
    double ycrossd;
    double zcrossd;
    double xcrossa;
    double ycrossa;
    double zcrossa;
    double dcons;
    double product;
    double torsionangle;
    //NCCcooridate[Naas][Natom][Nxyz]
    //Naas---i-1,i,i+1; Natom---N,Ca,C; Nxyz---x,y,z
    //NiCaiC'i---N=0, Ca=1; C=2;
    x10 = point1[0];
    y10 = point1[1];
    z10 = point1[2];
    x20 = point2[0];
    y20 = point2[1];
    z20 = point2[2];
    x30 = point3[0];
    y30 = point3[1];
    z30 = point3[2];
    cx1 = (y20-y10)*(z30-z10)-(y30-y10)*(z20-z10);
    cy1 = (x30-x10)*(z20-z10)-(x20-x10)*(z30-z10);
    cz1 = (x20-x10)*(y30-y10)-(x30-x10)*(y20-y10);
    //psi:NiCaiC'i--CaiC'iNi+1
    //phi:C'i-1NiCai--NiCaiC'i
    x1 = x20; 
    y1 = y20; 
    z1 = z20; 
    x2 = x30; 
    y2 = y30; 
    z2 = z30; 
    x3 = point4[0]; 
    y3 = point4[1]; 
    z3 = point4[2]; 
    cx2 = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
    cy2 = (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
    cz2 = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

    costheta = (cx1*cx2+cy1*cy2+cz1*cz2)/(sqrt(cx1*cx1+cy1*cy1+cz1*cz1)*sqrt(cx2*cx2+cy2*cy2+cz2*cz2));
    if ( costheta > 1 )
    {
        costheta = 1;
    }
    else if ( costheta < -1 )
    {
        costheta = -1;
    }
    torsionangle = acos(costheta);
    torsionangle = acos(costheta)*(180.0/3.14159);


    //direction
    nx = (x30-x20);
    ny = (y30-y20);
    nz = (z30-z20);
    dcons = nx*x20+ny*y20+nz*z20;

    //the cross point between line LD and plan BPBC: Dx(xcrossd,ycrossd,zcrossd)
    if ( ny != 0 )
    {
        ycrossd = dcons - nx*x3 - nz*z3 + y3/ny*(nx*nx+nz*nz);
        ycrossd = ycrossd/(ny+(nx*nx+nz*nz)/ny);
        xcrossd = (nx/ny)*(ycrossd-y3)+x3;
        zcrossd = (nz/ny)*(ycrossd-y3)+z3;
    }
    else if (nx != 0)
    {
        xcrossd = dcons - ny*y3 - nz*z3 + x3/nx*(ny*ny+nz*nz);
        xcrossd = xcrossd/(nx+(ny*ny+nz*nz)/nx);
        ycrossd = (ny/nx)*(xcrossd-x3)+y3;
        zcrossd = (nz/nx)*(xcrossd-x3)+z3;
    }
    else if (nz != 0)
    {
        zcrossd = dcons - nx*x3 - ny*y3 + z3/nz*(nx*nx+ny*ny);
        zcrossd = zcrossd/(nz+(ny*ny+nx*nx)/nz);
        ycrossd = (ny/nz)*(zcrossd-z3)+y3;
        xcrossd = (nx/nz)*(zcrossd-z3)+x3;
    }
	else /*(nx == 0 && ny == 0 && nz == 0)*/
	{
		/*avoid uninitialized value 2008-04-23, Nanjiang*/
		fprintf(stderr,"Error! Points 2 and 3 are overlapping!\n");
		return -999999;
	}

    //the cross point between line LA and plan BPBC: Ax(xcrossa,ycrossa,zcrossa)
    if ( ny != 0 )
    { 
        ycrossa = dcons - nx*x10 - nz*z10 + y10/ny*(nx*nx+nz*nz);
        ycrossa = ycrossa/(ny+(nx*nx+nz*nz)/ny);
        xcrossa = (nx/ny)*(ycrossa-y10)+x10;
        zcrossa = (nz/ny)*(ycrossa-y10)+z10;
    } 
    else if (nx != 0)
    { 
        xcrossa = dcons - ny*y10 - nz*z10 + x10/nx*(ny*ny+nz*nz);
        xcrossa = xcrossa/(nx+(ny*ny+nz*nz)/nx);
        ycrossa = (ny/nx)*(xcrossa-x10)+y10;
        zcrossa = (nz/nx)*(xcrossa-x10)+z10;
    }
    else if (nz != 0)
    {
        zcrossa = dcons - nx*x10 - ny*y10 + z10/nz*(nx*nx+ny*ny);
        zcrossa = zcrossa/(nz+(ny*ny+nx*nx)/nz);
        ycrossa = (ny/nz)*(zcrossa-z10)+y10;
        xcrossa = (nx/nz)*(zcrossa-z10)+x10;
    }
	else /*(nx == 0 && ny == 0 && nz == 0)*/
	{
		/*avoid uninitialized value 2008-04-23, Nanjiang*/
		fprintf(stderr,"Error! Points 2 and 3 are overlapping!\n");
		return -999999;
	}

    //The direction NABDx = AxB x BDx
    nx1 = (y20-ycrossa)*(zcrossd-z20)-(z20-zcrossa)*(ycrossd-y20);
    ny1 = (z20-zcrossa)*(xcrossd-x20)-(x20-xcrossa)*(zcrossd-z20);
    nz1 = (x20-xcrossa)*(ycrossd-y20)-(y20-ycrossa)*(xcrossd-x20);
    product = (nx*nx1+ny*ny1+nz*nz1)/sqrt(nx*nx+ny*ny+nz*nz)/sqrt(nx1*nx1+ny1*ny1+nz1*nz1);
    if (product<-0.99)
    { }
    else if (product>0.9999)
    { torsionangle=-1*torsionangle; }
    else
    {
    }
    return torsionangle;
}
/*}}}*/
float readfloatfromfile(char *bak,int k1,int k2)/*{{{*/
{
    int i;
    float number;
    char tbp[20];
    for (i=k1;i<=k2;i++)
    { tbp[i-k1] = bak[i];}
    tbp[i-k1] = '\0';
    number = float ( atof(tbp) );
    return number;
}
/*}}}*/
int readintfromfile(char *bak,int k1,int k2)/*{{{*/
{
    int i,number;
    char tbp[20];
    for (i=k1;i<=k2;i++)
    { tbp[i-k1] = bak[i];}
    tbp[i-k1] = '\0';
    number = atoi(tbp);
    return number;
}
/*}}}*/
int Aligment_Two_sequence_500(int *Seq1, int *Seq2, int Len1, int Len2, char *cseq1,char *cseq2, int *Align_Table,int *Max_Score)/*{{{*/
{
    int ik,ij,ip, Max_row, Max_col, Max_col_Length,Max_row_Length,Max_Length, Path_Vertical,Path_Horizontal, Path_Diagonal;
	int Max_score_last = -999999; /*initialize the Max_score_last, Nanjiang 2008-04-23*/

    int Result_Length, *Result_path,*Vertical_path;//Result_path for horizontal score,Vertical_path for vertical score of last column
    bool Gap_Prev_horizontal,*Gap_Prev_Row_Vertical;
    int Gap1_vertical,Gap2_horizontal, Gap_open, Gap1,Gap2,Score_Prev, Score_Upp_Left;
    int  NCode;
    char aas_str[] = "AVLIPFMKRHGSTCYNEWDQU  ";
    int Number_Run,loop_run,Each_width_len1, Max_row_part,NCount_eachpart;

    int Last_time_ij = 0;  /*do know which initial value should be set, ?????????, Nanjiang, 2008-04-23*/
	int Final_Length,N_tmp;
	int n_size;

    char *Result1_bak,*Result2_bak;
    int *Loop_beg ; /*2007-11-07, Nanjiang, changed to allocated memory*/ 
    int *Loop_end ;

#ifdef DEBUG
    Array1D <char> seq1char_1darray(Len1+1);
    Array1D <char> seq2char_1darray(Len2+2);
    char *seq1char = seq1char_1darray.array1D;
    char *seq2char = seq2char_1darray.array1D;
    for(ik = 0 ; ik < Len1; ik ++)
    {
        seq1char[ik] = aas_str[Seq1[ik]];
    }
    seq1char[ik] = '\0';

    for(ik = 0 ; ik < Len2; ik ++)
    {
        seq2char[ik] = aas_str[Seq2[ik]];
    }
    seq2char[ik] = '\0';
#endif

    strcpy(cseq1,"");/*added 2007-11-06, Nanjiang*/
    strcpy(cseq2,"");

    Gap_open = Align_Table[400];
    Gap1_vertical = Align_Table[401];
    Gap2_horizontal = Align_Table[402];

    Final_Length = 0;
    //check if the total space >=300*300; How many times should be run for each 300*300
    Number_Run = (Len1*Len2-1)/(500*500)+1;
    Each_width_len1 = Len1/Number_Run + 1;
    Array1D <int> Loop_beg_1darray(Number_Run+1);
    Array1D <int> Loop_end_1darray(Number_Run+1);
    Loop_beg = Loop_beg_1darray.array1D;
    Loop_end = Loop_end_1darray.array1D;

    for (loop_run=0; loop_run<Number_Run; loop_run++ )
    {
        Loop_beg[loop_run] = loop_run*Each_width_len1;
        Loop_end[loop_run] = (loop_run+1)*Each_width_len1+2;//one line is for cross-line for next loop_run, one line more just for "for (ik=1; ik<Loop_end[loop_run]; ik++)"
        if (   (Loop_end[loop_run]>=Len1) && (loop_run<(Number_Run-1))   )
        {
            Loop_end[loop_run] = Len1+1;
            Number_Run = loop_run + 1;
            break;
        }
        if (  loop_run == (Number_Run-1)  )
        {
            Loop_end[loop_run] = Len1+1;
        }
    }
    //set matrix
    Max_row = Len1 + 1;
    Max_col = Len2 + 1;
    Max_row_part = Each_width_len1 + 4;
    char  *Matrix_Aligment;
    n_size = Max_row_part*Max_col+1;
    Matrix_Aligment = new  char [n_size];
    for (ik=0; ik <n_size; ik++)
    {
        Matrix_Aligment[ik] = 0;
    }

    Max_row_Length = Max_row+10;
    Max_col_Length = Max_col+10;
    Max_Length = max(Max_row_Length,Max_col_Length)*2 + 50;
    Result_path = new int[Max_Length];
    Vertical_path = new int[Max_row_Length];
    Gap_Prev_Row_Vertical = new bool[Max_col_Length];
    Result1_bak = new char[Max_Length];
    Result2_bak = new char[Max_Length];


    for (loop_run=Number_Run-1; loop_run>=0; loop_run-- )//loop_run=0 is last part, loop_run<Number_Run-1 is the first part
    {

        NCount_eachpart = -1;
        Result_path[0] = 0;//here Result_path is used for storing scores
        Gap_Prev_Row_Vertical[0] = 1;
        Matrix_Aligment[0] = 0;
        //initialization of first row
        for (ij=1; ij<Max_col; ij++)
        {
            Result_path[ij] = Result_path[ij-1] + Gap2_horizontal;//horizontal
            Gap_Prev_Row_Vertical[ij] = 1;//first row
            Matrix_Aligment[ij] = 2;//horizontal
        }
        Vertical_path[0] = Result_path[Max_col-1];

        for (ik=1; ik<Loop_end[loop_run]; ik++)
        {
            if (   ik >= Loop_beg[loop_run]  )
            {
                NCount_eachpart++;
                if  (   (loop_run==0) && (NCount_eachpart==0)   )
                {
                    NCount_eachpart++;
                }
            }
            Score_Upp_Left = Result_path[0];
            Score_Prev = Result_path[0] + Gap1_vertical;
            Result_path[0] = Score_Prev;
            Gap_Prev_horizontal = 1;//first column

            if (  NCount_eachpart >= 0   )
            {
                Matrix_Aligment[NCount_eachpart*Max_col] = 0;//vertical
            }

            for (ij=1; ij<Max_col; ij++)
            {
                if (  Gap_Prev_horizontal == 0  )//gap open horizontal
                {
                    Gap2 = Gap_open;
                }
                else
                {
                    Gap2 = Gap2_horizontal;
                }
                if (  Gap_Prev_Row_Vertical[ij] == 0  )//gap open vertical
                {
                    Gap1 = Gap_open;
                }
                else
                {
                    Gap1 = Gap1_vertical;
                }

                Path_Horizontal = Score_Prev + Gap2;
                Path_Vertical = Result_path[ij] + Gap1;
                if (   (Seq1[ik-1]>=20) || (Seq2[ij-1]>=20)   )
                {
                    N_tmp = 403;
                }
                else
                {
                    N_tmp = Seq1[ik-1]*20+Seq2[ij-1];
                }
                Path_Diagonal = Score_Upp_Left + Align_Table[N_tmp];


                if (   (Path_Horizontal>=Path_Diagonal) && (Path_Horizontal>=Path_Vertical)   )
                {
                    Gap_Prev_horizontal = 1;
                    Gap_Prev_Row_Vertical[ij] = 0;
                    Max_score_last = Path_Horizontal;
                    NCode = 2;
                }
                else if (    (Path_Vertical>=Path_Diagonal) && (Path_Vertical>=Path_Horizontal)   )
                {
                    Gap_Prev_horizontal = 0;
                    Gap_Prev_Row_Vertical[ij] = 1;
                    Max_score_last = Path_Vertical;
                    NCode = 0;
                }
                else 
                {
                    Gap_Prev_horizontal = 0;
                    Gap_Prev_Row_Vertical[ij] = 0;
                    Max_score_last = Path_Diagonal;
                    NCode = 1;
                }


                Score_Upp_Left = Result_path[ij];
                Score_Prev = Max_score_last;
                Result_path[ij] = Max_score_last;

                if (  NCount_eachpart >= 0  )
                {
                    Matrix_Aligment[NCount_eachpart*Max_col+ij] = NCode;
                }

            }
            Vertical_path[ik] = Result_path[ij-1];

        }
        //Max_score_last = (*Matrix_Aligment)(Max_row,Max_col)

        Result_Length = 0;
        if (   loop_run == (Number_Run-1)  )//track back from the last one  Result_path[Max_col-1]
        {
            //trace back---calculate the "Result_path[]"; Path_Horizontal--0; Path_Vertical--1; Path_Diagonal--2
            ik = Loop_end[loop_run] - Loop_beg[loop_run];
            ij = Max_col;
            //first decide the last gaps
            if (  Max_score_last < Result_path[Max_col-2] )//gap for horizontal,,gap mark must be minus
            {
                ij = Max_col-2;
                while (   (Max_score_last<Result_path[ij]) && (ij>=1)   )
                {
                    Max_score_last = Result_path[ij];
                    Result_path[Result_Length] = 2;//horizontal gap
                    Result_Length++;
                    ij--;
                }
                ij = ij + 2;
            }
            else if (  Max_score_last < Vertical_path[Max_row-2] )//gap for vertical
            {
                ik = Loop_end[loop_run] - Loop_beg[loop_run] - 2;
                while (   (Max_score_last<Vertical_path[ik]) && (ik>=1)   )
                {
                    Max_score_last = Vertical_path[ik];
                    Result_path[Result_Length] = 0;
                    Result_Length++;
                    ik--;
                }
                ik = ik + 2;
            }
            Max_Score[0] = Max_score_last;

        }
        else//it is not needed for loop_run < (Number_Run-1) because of the position of Max_score is already set
        {
            ik = Loop_end[loop_run] - Loop_beg[loop_run] -1;//each time the ik is calculated one more time
            ij = Last_time_ij +1 ;//because after "while  (   (ik>=1) && (ij>=1)   )", there is --
            if (  ij > Max_col  )
            {
                ij = Max_col;
            }
        }

        //(*Matrix_Aligment)(ik,ij)--ik,ij begin at 1
        while  (   (ik>=1) && (ij>=1)   )
        {
            if (    Matrix_Aligment[(ik-1)*Max_col+ij-1] == 1     )//diagonal
            {
                Result_path[Result_Length] = 1;
                ik--;
                ij--;
            }
            else if (    Matrix_Aligment[(ik-1)*Max_col+ij-1] == 2     )//horizontal
            {
                Result_path[Result_Length] = 2;
                ij--;
            }
            else//vertical
            {
                Result_path[Result_Length] = 0;
                ik--;
            }
            Result_Length++;
        }
        Last_time_ij = ij;

        //check last
        if (  loop_run == 0  )//the ik and ij couldn't be zero at last loop
        {
            Result_Length--;//only one of ik and ij be zero
            if (  (ik>=1) && (ij==0)  )
            {

                for (ip=ik; ip>=1; ip--)
                {
                    Result_path[Result_Length] = Matrix_Aligment[(ip-1)*Max_col];
                    Result_Length++;
                }

            }
            else if (  (ij>1) && (ik==0)  )
            {
                for (ip=ij; ip>=1; ip--)
                {
                    Result_path[Result_Length] = Matrix_Aligment[ip-1];
                    Result_Length++;
                }

            }
            Last_time_ij = 0;
        }

        //set the aligment for Seq1 and Seq2, put result into Result1 and Result2
        Result_Length = Result_Length - 2;//(*Matrix_Aligment)(1,1)--- start point is taken away

        ik = Loop_beg[loop_run];//verical
        ij = Last_time_ij;//horinzontal
        if (  loop_run == 0  )
        {
            Result_Length++;//last
        }
        for (ip=Result_Length; ip>=0; ip--)
        {
            if (  Result_path[ip] == 1  )//diagonal
            {
                if (  (ik<Len1) && (ij<Len2)   )
                {
                    Result1_bak[Result_Length-ip] = aas_str[Seq1[ik]];
                    Result2_bak[Result_Length-ip] = aas_str[Seq2[ij]];
                    ik++; 
                    ij++;
                }
                else
                {
                }
            }
            else  if (  Result_path[ip] == 2  )//horinzontal
            {
                if (  ij < Len2   )
                {
                    Result1_bak[Result_Length-ip] = aas_str[21];//gap
                    Result2_bak[Result_Length-ip] = aas_str[Seq2[ij]];
                    ij++;
                }
                else
                {
                }
            }
            else//vertical
            {
                if (  ik < Len1   )
                {
                    Result1_bak[Result_Length-ip] = aas_str[Seq1[ik]];
                    Result2_bak[Result_Length-ip] = aas_str[21];//gap
                    ik++;
                }
                else
                {
                }
            }
        }

        Result1_bak[Result_Length+1] = '\0'; /*bug fixed, nanjiang, 2007-11-06, Result1_bak should be added to the front of cseq1*/
        Result2_bak[Result_Length+1] = '\0';
        strcat(Result1_bak, cseq1);
        strcat(Result2_bak, cseq2);
        strcpy(cseq1, Result1_bak);
        strcpy(cseq2, Result2_bak);


        //for (ik=0; ik<Result_Length+1; ik++)
        //{
        //    cseq1[Final_Length+ik] = Result1_bak[ik];[>bug fixed, nanjiang, 2007-11-06, Result1_bak should be added to the front of cseq1<]
        //    cseq2[Final_Length+ik] = Result2_bak[ik];
        //}
        Final_Length = Final_Length + Result_Length+1;
    } 
    if (  Final_Length >= Max_Length  )
    {
        return 0;
    }

    delete [] Matrix_Aligment; 
    delete [] Result_path;
    delete [] Gap_Prev_Row_Vertical;
    delete [] Vertical_path;
    delete [] Result1_bak;
    delete [] Result2_bak;
    cseq1[Final_Length] = '\0';
    cseq2[Final_Length] = '\0';
    return Final_Length;
}
/*}}}*/
int GetShapeString(const char *pdbfile, char* chainIDList, char *title, int *Align_Table, int8 *** Shapestring, bool isPrintAASeq, const char* aaSeqFile, const char* dsspFile, FILE *fpout, FILE *fpAlign)/*{{{*/
/*****************************************************************************
 * using SeqMatchAlign_Protein, 2007-11-07 17:49:06 Wednesday Week 44
 ****************************************************************************/
{
    char Catom[6] = "";
    char Caas[5] = "";

    int ik,Ntemp;

    int     AASseriesnumber1;
    int     AASseriesnumber2;
    int     AASseriesnumber3;
    int     Sseriesnumber;
    int     NCAC[3][3][3];
    int     NAAS;
    int     NATOM;
    int     im;
    int     ij;
    double  coordinates[3][3][3];
    double  Xuccupncy;
    double  Xuccupncy_all[3];                              //position-atom-xyz
    //char    caasstring[]         = "AVLIPFMKRHGSTCYHEWDQ";
    char    cshapestring[]       = "SRUVKATG-";
    char    CatominATOM[4];
    int     Len_Seq;
    int     Len_Seq_atom = 0;
    int     Nnew_long;
    int    *Seq_SEQ;
    int    *Seq_ATOM;
    int    *Series_Atom;
    bool    Read_Atom;
    float  *vphi;
    float  *vpsi;

    char    *cseq1 = new char [MAXSeqLen*2+1];
    char    *cseq2 = new char [MAXSeqLen*2+1];
    char    *cshape= new char [MAXSeqLen*2+1];/*setting the maxseqlen to a larger value, 2007-11-06, Nanjiang*/

    for(ik = 0  ; ik < MAXSeqLen*2+1; ik ++)
    {
        cshape[ik] = '-'; /*2007-11-09, initialize cshape to '-'*/
    }

    Seq_ATOM = new int [MAXSeqLen];
    Seq_SEQ =  new int [MAXSeqLen]; 
    Series_Atom = new int [MAXSeqLen];
    vphi = new float[MAXSeqLen];
    vpsi = new float[MAXSeqLen];


    int  nbeg;
    int  nend;
    int  first_O = 0; /*don't know which initial value to be set, 2008-04-23, Nanjiang*/

    Array1D <char> aaSeq_1darray(MAXSeqLen+1);
    char *aaSeq = aaSeq_1darray.array1D;
    Array1D <char> atomSeq_1darray(MAXSeqLen+1);
    char *atomSeq = atomSeq_1darray.array1D;

    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char chainID=chainIDList[0];

    Len_Seq = GetSeq_SEQRES(pdbfile, chainID, aaSeq) ;

    if (isUsingPDBAtom)/*{{{*/
    {
        FILE *fpread = NULL;
        fpread = fopen(pdbfile,"r");
        if (  fpread == NULL  ) { return -1; }

        Nnew_long = 0;
        Len_Seq = 0;
        while((linesize  = fgetline(fpread, line, maxline)) != EOF)
        //while (  fgets(line,900,fpread) != NULL  )
        {
            if (   linesize < 20  )
            {
                continue;
            }
            for (ik=0; ik<4; ik++)
            {
                Catom[ik] = line[ik];
            }
            Catom[ik] = '\0';
            if (  (strcmp(Catom,"SEQR")==0) && (line[11]==chainID)  )  /*read in sequence from SEQRES record, 2007-10-21, nanjiang*//*{{{*/
            {
                if (  Nnew_long == 0  )
                {
                    Nnew_long = readintfromfile(line,13,16);
                    //Seq_ATOM = new int[Nnew_long+50];
                    //Seq_SEQ = new int[Nnew_long+50];
                    //Series_Atom = new int[Nnew_long+50];
                    //vphi = new float[Nnew_long+50];
                    //vpsi= new float[Nnew_long+50];
                    for (ik=0; ik<Nnew_long; ik++)
                    {
                        Seq_ATOM[ik] = 20;
                        Seq_SEQ[ik] = 20;
                        vphi[ik] = -360;
                        vpsi[ik] = -360;
                    }
                }

                for (ik=19; ik<71; ik=ik+4)   /*originally, ik < 81, bug fixed by nanjiang 2007-07-12*/
                {
                    for (ij=ik; ij<ik+3; ij++)
                    {
                        Caas[ij-ik] = line[ij];
                    }
                    Caas[ij-ik] = '\0';

                    if(strcmp(Caas, "   ") == 0)    /*to the end of sequence, neglect the trailing blanks
                                                      bug fixed 2007-10-21, the last blank Caas will be "   ", not ""*/
                    {
                        continue;
                    }

                    for (ij=0; ij<20; ij++)
                    {
                        if (  strcmp(Caas,AAs[ij]) == 0 )
                        {
                            break;
                        }
                    }
                    Seq_SEQ[Len_Seq] = ij;
                    Len_Seq++;
                }
            }/*}}}*/
            else if ((strcmp(Catom,"ATOM")==0) && (line[21]==chainID)  )/*{{{*/  /*read in the sequence form ATOM record, 2007-10-21, nanjiang*/
            {
                //******************************
                AASseriesnumber1 = -100;
                AASseriesnumber2 = -100;
                AASseriesnumber3 = -100;
                for (ik=0; ik<3; ik++)
                {
                    for (ij=0; ij<3; ij++)
                    {
                        for (im=0; im<3; im++)
                        {
                            coordinates[ik][ij][im] = 0;
                            NCAC[ik][ij][im] = 0;
                        }
                    }
                    Xuccupncy_all[ik] = 0;
                }
                Read_Atom = 1;
                Len_Seq_atom = 0;
                while (  Read_Atom  )//(  (strcmp(Catom,"ATOM")==0) && (line[21]==chainID)  )
                {
                    for (ik=0; ik<4; ik++)
                    {
                        Catom[ik] = line[ik];
                    }
                    Catom[ik] = '\0';
                    if (  (line[21]!=chainID) || (line[0]=='T')  )
                    {
                        break;
                    }
                    else if (  (line[21]==chainID) && (strcmp(Catom,"ATOM")!=0)   )//ANISOU and HETATM
                    {
                        if ( fgets(line,900,fpread) == NULL )
                        {
                            break ;
                        }
                        continue;
                    }

                    //read series number
                    Sseriesnumber = readintfromfile(line,22,25);
                    //read AAS name
                    for (ik=0; ik<3; ik++)
                    {
                        Caas[ik] = line[17+ik];
                    }
                    Caas[ik] = '\0';
                    for (ik=0; ik<20; ik++)
                    {
                        if (  strcmp(Caas,AAs[ik]) == 0 )
                        {
                            break;
                        }
                    }
                    NAAS = ik;
                    if (  Sseriesnumber != AASseriesnumber3  )
                    {
                        Seq_ATOM[Len_Seq_atom] = NAAS;
                        Series_Atom[Len_Seq_atom] = Sseriesnumber;
                        Len_Seq_atom++;
                        AASseriesnumber3 = Sseriesnumber;
                    }
                    //read atom name
                    for (ik=0; ik<4; ik++)
                    {
                        CatominATOM[ik] = line[12+ik];
                    }
                    Ntemp = 0;
                    for (ik=0; ik<4; ik++)
                    {
                        if (  CatominATOM[ik] != ' ' )
                        {
                            CatominATOM[Ntemp] = CatominATOM[ik];
                            Ntemp++;
                        }
                    }
                    if (  Ntemp == 1 )
                    {
                        CatominATOM[Ntemp] = '\0';
                    }
                    else
                    {
                        CatominATOM[2] = '\0';
                    }

                    if (  strcmp(CatominATOM,"N") == 0  )
                    {
                        NATOM = 0;
                    }
                    else if (  strcmp(CatominATOM,"CA") == 0  )
                    {
                        NATOM = 1;
                    }
                    else if (  strcmp(CatominATOM,"C") == 0  )
                    {
                        NATOM = 2;
                    }
                    else if (  strcmp(CatominATOM,"O") == 0  )
                    {
                        NATOM = 3;
                    }
                    else
                    {
                        NATOM = 5;
                    }
                    if (  NATOM == 5  )
                    {
                        if ( fgets(line,900,fpread) == NULL )
                        {
                            break ;
                        }
                        if (  line[21] != chainID  )
                        {
                            Read_Atom = 0;
                        }
                        continue;
                    }
                    if (  NATOM <= 2 )//N, CA, C
                    {
                        //occupnacy
                        Xuccupncy = readfloatfromfile(line,54,59);//occupnacy
                        coordinates[2][NATOM][0] = coordinates[2][NATOM][0] + Xuccupncy*readfloatfromfile(line,30,37);//x
                        coordinates[2][NATOM][1] = coordinates[2][NATOM][1] + Xuccupncy*readfloatfromfile(line,38,45);//y
                        coordinates[2][NATOM][2] = coordinates[2][NATOM][2] + Xuccupncy*readfloatfromfile(line,46,53);//z
                        Xuccupncy_all[NATOM] =  Xuccupncy_all[NATOM] + Xuccupncy;
                        first_O = 1;
                    }
                    else if ( first_O == 1 )//O and first_O == 1
                    {
                        //if (  ((AASseriesnumber3-AASseriesnumber2)==1) && ((AASseriesnumber2-AASseriesnumber1)==1)  )
                        if (  Len_Seq_atom >= 3     )
                        {
                            for (ik=0; ik<3; ik++)
                            {
                                for (ij=0; ij<3; ij++)
                                {
                                    if (  NCAC[ik][ij] == 0  )
                                    {
                                        ik = 10; 
                                        break;
                                    }
                                }
                            }
                            for (ij=0; ij<3; ij++)
                            {
                                if (  Xuccupncy_all[ij] < 0.1 )
                                {
                                    ik = 10; 
                                    break;
                                }
                            }
                            if (  ik < 10 )//torsion angles
                            {
                                //recalculate the last residue base on the Xuccupncy_all[3]
                                for (ij=0; ij<3; ij++)
                                {
                                    for (im=0; im<3; im++)
                                    {
                                        coordinates[2][ij][im] = coordinates[2][ij][im]/Xuccupncy_all[ij];
                                    }
                                }
                                //torsion angles: 
                                //phi--C'i-1NiCai--NiCaiC'i,   C'i-1NiCaiC'i---1,2,3,4
                                vphi[Len_Seq_atom-2] = float ( Calculate_torsion_angles_phi_psi(coordinates[0][2],coordinates[1][0],coordinates[1][1],coordinates[1][2]) );
                                //phi--NiCaiC'i--CaiC'iNi+1,   NiCaiC'iNi+1---1,2,3,4
                                vpsi[Len_Seq_atom-2] = float ( Calculate_torsion_angles_phi_psi(coordinates[1][0],coordinates[1][1],coordinates[1][2],coordinates[2][0]) );
                            }
                        }
                        AASseriesnumber1 = AASseriesnumber2;
                        AASseriesnumber2 = AASseriesnumber3;
                        for (ik=0; ik<2; ik++)
                        {
                            for (ij=0; ij<3; ij++)
                            {
                                for (im=0; im<3; im++)
                                {
                                    coordinates[ik][ij][im] = coordinates[ik+1][ij][im];
                                    NCAC[ik][ij][im] = NCAC[ik+1][ij][im];
                                }
                            }
                        }
                        for (ij=0; ij<3; ij++)
                        {
                            for (im=0; im<3; im++)
                            {
                                NCAC[2][ij][im] = 0;
                                coordinates[2][ij][im] = 0;
                            }
                            Xuccupncy_all[ij] = 0;
                        }
                        first_O = 0;
                    }//O atom
                    if ( fgets(line,900,fpread) == NULL )
                    {
                        break ;
                    }
                    if (  line[21] != chainID  )
                    {
                        Read_Atom = 0;
                    }
                }//while

                //********************************
            }//if/*}}}*/

        }
        fclose(fpread);

        for(ik = 0 ; ik < Len_Seq_atom; ik++)
        {
            atomSeq[ik] = AA3To1(AAs[Seq_ATOM[ik]] );
        }
        atomSeq[ik] = '\0';
        for( ik = 0 ; ik < Len_Seq; ik ++)
        {
            aaSeq[ik] = AA3To1(AAs[Seq_SEQ[ik]] );
        }
        aaSeq[ik] = '\0';

        //shape line
        for (ik=0; ik<Len_Seq_atom; ik++)
        {
            cshape[ik] = '-';
            if (    (vphi[ik]>=-180) && (vphi[ik]<=180) && (vpsi[ik]>=-180) && (vpsi[ik]<=180)    )
            {
                nbeg = int ( (vphi[ik]+180)/2 );
                nend = int ( (vpsi[ik]+180)/2 );
                if (  Seq_ATOM[ik] < 20  )
                {
                    cshape[ik] = cshapestring[Shapestring[Seq_ATOM[ik]][nbeg][nend]];
                }
            }
        }
        cshape[ik] = '\0';
    }/*}}}*/
    else /*getting torsion angles from DSSP file*//*{{{*/
    {
        char pdbid[SIZE_PDBID+1] = "";
        char id[SIZE_CHAIN_ID+1] = "";
        char tmpdsspfile[MAX_PATH+1] = "";
        int sizePDBFile = strlen(pdbfile);
        my_strcpy(pdbid, pdbfile+sizePDBFile-8, SIZE_PDBID);
        my_strupr(pdbid);
        if (strcmp(dsspFile, "") == 0)
        {
            if (GetDSSPFilePath(pdbid, tmpdsspfile) == NULL)
            {
                return -1; /*invoking dssp file has not been implemented*/
            }
        }
        else
        {
            my_strcpy(tmpdsspfile, dsspFile, MAX_PATH);
#ifdef DEBUG
            printf("dsspFile=%s\n", dsspFile);
#endif
        }
        
        sprintf(id,"%s%c", pdbid, chainID);
        
        Chain dsspChain;
        InitChain(&dsspChain);
        dsspChain.aaSeq = new char[MAXSeqLen+1] ;
        dsspChain.phi = new float[MAXSeqLen+1] ;
        dsspChain.psi = new float[MAXSeqLen+1] ;
        if (GetDSSPChain(id, &dsspChain, tmpdsspfile) < 0)
        {
            return -1;
        }
        Len_Seq_atom = dsspChain.numRes;
#ifdef DEBUG
            printf("dsspChain.numRes=%d\n", dsspChain.numRes);
            for(ik=0;ik<dsspChain.numRes;ik++)
            {
                printf("psi[%d]=%f \t phi[%d]=%f\n", ik, dsspChain.psi[ik], ik, dsspChain.phi[ik]);
            }
#endif
        for(ik = 0 ; ik < Len_Seq_atom; ik ++)
        {
            if(islower(dsspChain.aaSeq[ik]))
            {
                dsspChain.aaSeq[ik] = 'C';
            }
        }
        my_strcpy(atomSeq, dsspChain.aaSeq, MAXSeqLen);

        for (ik=0; ik<Len_Seq_atom; ik++)
        {
            for (ij=0; ij<20; ij++)
            {
                if ( atomSeq[ik] == AAalphabet[ij] ) { break; }
            }
            Seq_ATOM[ik] = ij;
        }
        
        for(ik = 0 ; ik < Len_Seq_atom; ik ++)
        {
            vphi[ik] = dsspChain.phi[ik];
            vpsi[ik] = dsspChain.psi[ik];
        }
        
        DeleteChain(&dsspChain);
    }/*}}}*/

#ifdef DEBUG
    fprintf(stdout, ">%s length = %d sequence from SEQRES record\n", title, Len_Seq);
    char seq_seqres[5000] = "";
    int i= 0;
    for( i = 0 ; i < Len_Seq; i ++)
    {
        seq_seqres[i] = AA3To1(AAs[Seq_SEQ[i]] );
    }
    seq_seqres[i] = '\0';
    WriteFastaSeq(seq_seqres, stdout);
#endif


    if(strcmp(aaSeqFile,"") != 0) /*if pdbaa file supplied*//*{{{*/
    {
        int seqtype = 0;
        if ((Len_Seq = ReadSeq_FASTA(aaSeqFile, aaSeq,&seqtype, MAXSeqLen ))<= 0)
        {
            fprintf(stderr,"Read sequence file %s failed!\n",aaSeqFile);
            return -1;
        }
        else
        {
            for(ik = 0 ; ik < Len_Seq; ik++)
            {
                Seq_SEQ[ik] = Char2Digit(aaSeq[ik], AAalphabet, sizeAAalpabet);
            }
        }
    }/*}}}*/

#ifdef DEBUG
    char aaSeqTMP1[4000];
    int tmpii;
    for(tmpii = 0 ; tmpii < Len_Seq; tmpii++)
    {
        aaSeqTMP1[tmpii] = AA3To1(AAs[Seq_SEQ[tmpii]] );
    }
    aaSeqTMP1[tmpii] = '\0';
    fprintf(stdout,">Seq_SEQ\n");
    WriteFastaSeq(aaSeqTMP1, stdout);

    for(tmpii = 0 ; tmpii < Len_Seq_atom; tmpii++)
    {
        aaSeqTMP1[tmpii] = AA3To1(AAs[Seq_ATOM[tmpii]] );
    }
    aaSeqTMP1[tmpii] = '\0';

    fprintf(stdout,">Seq_ATOM\n");
    WriteFastaSeq(aaSeqTMP1, stdout);
#endif

    /*this code is added by nanjiang, 2007-11-07 */
    Array1D <char> shString_seq_1darray(Len_Seq+1);
    char *shString_seq = shString_seq_1darray.array1D;
    char title1[20] = "seq";
    char title2[20] = "atomseq";
    int length1 = Len_Seq;
    int length2 = Len_Seq_atom;
    Array1D <char>  alignXstr_1darray(length1 + length2 + 2); 
    Array1D <char>  alignYstr_1darray(length1 + length2 + 2);
    Array1D <int>   alignRel_1darray (length1 + length2 + 1);
    char *alignXstr = alignXstr_1darray.array1D;
    char *alignYstr = alignYstr_1darray.array1D;
    int  *alignRel  = alignRel_1darray.array1D;
    int alignLength = 0;
    bool isPrintToScreen = false;
    if (fpAlign != NULL)
    {
        fprintf(fpAlign,"***************************************************\n\n");
        fprintf(fpAlign,"alignment between aaseq and atomSeq for id : %s\n",title);
    }
    //if (isPrintToScreen)
    //{
    //    fprintf(stdout,"***************************************************\n\n");
    //    fprintf(stdout,"alignment between whole chain seq and scop seq for id : %s\n",title);
    //}
    alignLength = SeqMatchAlign_Protein( aaSeq, atomSeq,title1,title2, alignXstr,alignYstr,alignRel,fpAlign, isPrintToScreen);
    fflush(fpAlign);
    int seq_chain_cnt= 0 ;
    int atom_seq_cnt = 0 ;
    int cntIDT = 0;
    int j;
    Array1D <float> seqPhi_1darray(Len_Seq);
    Array1D <float> seqPsi_1darray(Len_Seq);
    float *seqPhi = seqPhi_1darray.array1D;
    float *seqPsi = seqPsi_1darray.array1D;
    /* post alignment processing*/
    for(j = 0 ; j < alignLength ; j ++)/*{{{*/
    {
        if(alignRel[j] == IDT)
        {
            //shString_seq[seq_chain_cnt] = cshape[atom_seq_cnt];
            seqPhi[seq_chain_cnt] = vphi[atom_seq_cnt];
            seqPsi[seq_chain_cnt] = vpsi[atom_seq_cnt];
            seq_chain_cnt ++;
            atom_seq_cnt ++;
            cntIDT ++;
        }
        else if(alignRel[j] == GAP)
        {
            if(alignYstr[j] == CHAR_INDEL)
            {
                //shString_seq[seq_chain_cnt] = '-' ; //set resSer to INIT_RESSEQ, if this residue does not exist in ATOM record
                seqPhi[seq_chain_cnt] = INIT_FLOAT;
                seqPsi[seq_chain_cnt] = INIT_FLOAT;
                seq_chain_cnt ++;
            }
            else if(alignXstr[j] == CHAR_INDEL) // if there are more residues in ATOM record than in SEQRES record, e.g. in 1EJG
            {
                atom_seq_cnt ++;
                //	fprintf(fpLog,"false align, the SEQRES missing residue\n");
                //	exit(0);
            }
        }
        else if(alignRel[j] == MIS)
        {
            //shString_seq[seq_chain_cnt] = cshape[atom_seq_cnt];
            seqPhi[seq_chain_cnt] = vphi[atom_seq_cnt];
            seqPsi[seq_chain_cnt] = vpsi[atom_seq_cnt];
            seq_chain_cnt ++;
            atom_seq_cnt ++;
        }
    }/*}}}*/

    /*calculate shape string for the aaSeq, from the phi-psi angles */
    for (ik=0; ik<Len_Seq; ik++)
    {
        shString_seq[ik] = '-';
        if (    (seqPhi[ik]>=-180.0) && (seqPhi[ik]<=180.0) && (seqPsi[ik]>=-180.0) && (seqPsi[ik]<=180.0)    )
        {
            nbeg = int ( (seqPhi[ik]+180.0)/2 );
            nend = int ( (seqPsi[ik]+180.0)/2 );
            int daa = 0;
            for (ij=0; ij<20; ij++) { if (aaSeq[ik] == AAalphabet[ij] ) { break; } }
            daa = ij;

            if (  daa < 20  )
            {
                shString_seq[ik] = cshapestring[Shapestring[daa][nbeg][nend]];
            }
        }
    }
    shString_seq[Len_Seq] = '\0';


    /*write out the shString_seq, or sequence*/
    fprintf(fpout,">%s  %d\n",title,Len_Seq);/*set the length to the sequence length, Nanjiang, 2007-11-06 */


    if(!isPrintAASeq)
    {
        WriteFastaSeq(shString_seq, fpout, 0, Len_Seq, printLineLength);
    }
    else 
    {
        Array1D <char> tmpstr_1darray(printLineLength+1);
        char *tmpstr = tmpstr_1darray.array1D;
        int beg = 0;
        while(beg <  Len_Seq)
        {
            my_strcpy(tmpstr, aaSeq + beg, printLineLength);
            fprintf(fpout, "%s\n", tmpstr);
            my_strcpy(tmpstr, shString_seq + beg, printLineLength);
            fprintf(fpout, "%s\n", tmpstr);
            beg += printLineLength;
            fprintf(fpout, "\n" );
        }
    }

    delete [] cseq1  ;
    delete [] cseq2  ;
    delete [] cshape ;
    delete [] Seq_ATOM;
    delete [] Seq_SEQ;
    delete [] Series_Atom;
    delete [] vphi;
    delete [] vpsi;

    return Len_Seq;
} /*}}}*/

int main( int argc, char *argv[])/*{{{*/
{
    if ( argc < 2 )
    {
        PrintHelp();
        return -1;
    }
    char outfile[MAX_PATH +1 ] = "";
    char alignFile[MAX_PATH +1 ] = "";
    char listfile[MAX_PATH +1 ] = "";
    bool isPrintAASeq = false;
    char pdbfile[MAX_PATH+1] = "";
    char aaSeqFile[MAX_PATH+1] = "";
    char dsspFile[MAX_PATH+1] = "";
    char chainIDList[300] = ""; /*changed 2009-12-21*/
    char title[300] = "";
    int max_num_tmpstr = 10;
    int max_size_tmpstr = 200;
    Array2D <char> tmpstrs_2darray(max_num_tmpstr, max_size_tmpstr+1);
    char ** tmpstrs= tmpstrs_2darray.array2D;
    bool isReadTextRawdata = false;

    char datadir[MAX_PATH+1]       = "";
    GetDataDir(datadir);
    char shapestring_rawdata_path[MAX_PATH+1] = "";
    sprintf(shapestring_rawdata_path,"%s/%s/%s",datadir,"shapestring", "rawdata");


    int i = 1;
    int j = 0;

    while( i < argc)/*{{{*/
    {
        if(strcmp(argv[i] , "-h") == 0 ||strcmp(argv[i], "--help") == 0)
        {
            PrintHelp();
            return 0;
        }
        else if (strcmp(argv[i], "-l") == 0)
        {
            my_strcpy(listfile, argv[i+1], MAX_PATH);
            i += 2;
        }
        else if (strcmp(argv[i], "-o") == 0)
        {
            my_strcpy(outfile, argv[i+1], MAX_PATH);
            i += 2;
        }
        else if (strcmp(argv[i], "--align") == 0)
        {
            my_strcpy(alignFile, argv[i+1], MAX_PATH);
            i += 2;
        }
        else if (strcmp(argv[i], "-a") == 0)
        {
            isPrintAASeq = true;
            i ++;
        }
        else if (strcmp(argv[i], "--nocase") == 0)
        {
            isChainIDCaseSensitive = false;
            i ++;
        }
        else if (strcmp(argv[i], "--rtxt") == 0)
        {
            isReadTextRawdata = true;
            i ++;
        }
        else if (strcmp(argv[i], "--max-seq") == 0)
        {
            MAXSeqLen = atoi(argv[i+1]);
            i += 2;
        }
        else if (strcmp(argv[i], "--forcedssp") == 0)
        {
            isForceRunDSSP = true; /*force running dssp program to generate dssp file*/
            i ++;
        }
        else if (strcmp(argv[i], "--pdb") == 0)
        {
            isUsingPDBAtom = true; /*using PDB Atom coordinates to calculate torsion angles*/
            i ++;
        }
        else if (strcmp(argv[i], "--linesize") == 0)
        {
            printLineLength = atoi(argv[i+1]);
            i += 2;
        }
        else if (strcmp(argv[i], "--data") == 0)
        {
            my_strcpy(shapestring_rawdata_path, argv[i+1], MAX_PATH);
            i += 2;
        }
        else
        {
            if(j < max_num_tmpstr)
            {
                my_strcpy(tmpstrs[j], argv[i], max_size_tmpstr);
                j ++;
            }
            i ++;
        }
    }/*}}}*/
    int numPara = j;



    if( numPara >= 1)
    {
        my_strcpy(pdbfile, tmpstrs[0], MAX_PATH);
    }
    if(numPara >= 2)
    {
        if(tmpstrs[1][1] != '\0')
        {
            my_strcpy(chainIDList, &(tmpstrs[1][1]), 300) ; /*chainIDList == "", means print all chains*/
            if (!isChainIDCaseSensitive)
            {
                my_strupr(chainIDList);
            }
            //if (!isChainIDCaseSensitive){chainID = toupper(tmpstrs[1][1]);}
            /*2008-07-07, set an option to control the case sensitivity of the
             * chainID*/
            //else { chainID = tmpstrs[1][1]; } [>2008-09-07, when isChainIDCaseSensitive is true, chainID should be assigned as well<]
        }
    }
    if (numPara >= 3)
    {
        my_strcpy(title, tmpstrs[2], 200);
    }
    if (numPara >= 4)
    {
        my_strcpy(aaSeqFile, tmpstrs[3], 200);
    }
    if (numPara >= 5)
    {
        my_strcpy(dsspFile, tmpstrs[4], 200);
    }

    if(strcmp(title, "") == 0)
    {
        my_strcpy(title, pdbfile, 200);
    }


    if(strcmp(listfile, "") == 0 && strcmp(pdbfile, "") == 0)
    {
        fprintf(stderr,"Error! Neither listfile nor pdbfile are set in the argument list\n");
        return -1;
    }

    FILE *fpout = NULL;

    if (  strcmp(outfile, "") != 0 )
    {
        fpout = fopen(outfile,"w");
        checkfilestream(fpout, outfile, "w");
    }
    else
    {
        fpout = stdout ;
    }

    FILE *fpAlign = NULL;  /*added 2007-11-07, Nanjiang*/
    if( strcmp(alignFile,"") != 0)
    {
        if(strcasecmp(alignFile, "stdout") == 0)
        {
            fpAlign = stdout;
        }
        else
        {
            fpAlign = fopen(alignFile, "w");
            checkfilestream(fpAlign,alignFile,"w");
        }
    }


    //get AAalphabet
    for (i = 0 ;i  < sizeAAalpabet; i ++)
    {
        AAalphabet[i] = AA3To1(AAs[i]);
    }
    AAalphabet[i] = '\0';


    Array1D <int> Align_Table_1darray(20*20+100);
    int *Align_Table = Align_Table_1darray.array1D;

    Array3D <int8> Shapestring_3darray (20,181,181);
    int8 *** Shapestring = Shapestring_3darray.array3D;

#ifdef DEBUG_TIME
    clock_t start, finish;
    double duration;
    int tmp;
    start = clock();
#endif

    GetParameter(shapestring_rawdata_path, Shapestring, Align_Table, isReadTextRawdata);

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"reading rawdata cost %lf seconds\n", duration);
#endif


#ifdef DEBUG_TIME
    start = clock();
#endif

    if (strcmp(pdbfile, "") != 0)
    {
        GetShapeString(pdbfile, chainIDList, title,  Align_Table, Shapestring, isPrintAASeq, aaSeqFile, dsspFile, fpout, fpAlign);
    }
    else // if(strcmp(listfile,"") != 0)
    {

        FILE *fplist = fopen(listfile, "r");
        checkfilestream(fplist, listfile, "r");
        int linesize;
        int maxline = 300;
        Array1D <char> line_1darray(maxline+1);
        char *line = line_1darray.array1D;
        char first_non_blank_char = ' ';

        while((linesize = fgetline(fplist, line, maxline)) != EOF)
        {
            first_non_blank_char = ' ';
            sscanf(line, " %c", &first_non_blank_char);
            if(linesize <= 0 || first_non_blank_char == '#') continue;
            char *linetrimed = strtrim(line);
            char *pch ;
            pch = strtok (linetrimed, WHITE_SPACE);
            int j = 0;
            strcpy(pdbfile, "");
            strcpy(chainIDList,"") ;
            strcpy(title, "");
            strcpy(aaSeqFile,"");
            strcpy(dsspFile,"");
            while(pch != NULL)
            {
                if(j == 0) { my_strcpy(pdbfile, pch, MAX_PATH); }
                else if( j == 1)
                {
                    if(*(pch+1) != '\0')
                    {
                        my_strcpy(chainIDList,  (pch+1), 300); 
                        if (!isChainIDCaseSensitive) { my_strupr(chainIDList); }
                    }
                        
                }
                else if(j == 2) { my_strcpy(title, pch, 200); }
                else if(j == 3) { my_strcpy(aaSeqFile, pch, MAX_PATH); }
                else if(j == 4) { my_strcpy(dsspFile, pch, MAX_PATH); }
                j ++;
                pch = strtok(NULL, WHITE_SPACE);
            }
            if(strcmp(title, "") == 0)
            {
                my_strcpy(title, pdbfile, 200);
            }

            GetShapeString(pdbfile, chainIDList, title,  Align_Table, Shapestring, isPrintAASeq, aaSeqFile,dsspFile,fpout, fpAlign);
        }
    }

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"creating shapestring cost %lf seconds\n", duration);
#endif

    if(fpout != NULL && fpout != stdout ) fclose(fpout);
    if(fpAlign != NULL && fpAlign != stdout ) fclose(fpAlign);
    return 0;
}

