#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "array.h"
#include "myfunc.h"

#undef MAX_NUM_PDBID
#define MAX_NUM_PDBID 90000

#undef SIZE_PDBID
#define SIZE_PDBID 4

void PrintHelp()
{
    fprintf(stdout,"Usage: getpdbfilepath [options] pdbid1 pdbid2 ...\n");
    fprintf(stdout,"  pdbid is case insensitive, first 4 characters are used for each argv\n");
    fprintf(stdout,"  piping supported" );
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  -d  <dir>  : path for main files\n"); 
    fprintf(stdout,"             : default = $DATADIR/pdb/data/structures/divided/pdb_dcp\n");
    fprintf(stdout,"  -do <dir>  : path for obsolete files\n");
    fprintf(stdout,"             : default = $DATADIR/pdb/data/structures/obsolete/pdb_dcp\n");
    fprintf(stdout,"  -md  <dir> : path for model structures\n");
    fprintf(stdout,"             : default = $DATADIR/pdb/data/structures/models/current/pdb_dcp\n");
    fprintf(stdout,"  -mdo <dir> : path for obsolete model structures\n");
    fprintf(stdout,"             : default = $DATADIR/pdb/data/structures/models/obsolete/pdb_dcp\n");
    fprintf(stdout,"  -l  <file> : pdbid file list\n");
    fprintf(stdout,"  -h|--help  : print this help message\n");

    fprintf(stdout,"Nanjiang Shu, created 2007-01-30, updated 2011-03-17\n ");
}

int ReadInPDBIDs(FILE *fp, char **pdbids, int start_idx)/*{{{*/
    //read in the chain ids from the file stream,
{
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int i = 0;
    while((linesize = fgetdelim(fp, line, WHITE_SPACE, maxline))!= EOF)
    {
        if(linesize > 0)
        {
            my_strcpy(pdbids[start_idx+i], line , SIZE_PDBID);
            i ++;
        }
    }
    return start_idx + i;
}
/*}}}*/
int main(int argc , char** argv)/*{{{*/
{
    //if ( argc < 2)
    //{
    //    PrintHelp();
    //    return 1;
    //}

    int i;
    char pdbpath[MAX_PATH+1] = "";
    char pdbobsoletepath[MAX_PATH+1] = "";
    char pdbmodelspath[MAX_PATH+1] = "";
    char pdbmodelsobsoletepath[MAX_PATH+1] = "";
    char pdbidlistfile[MAX_PATH+1] = "";
    char outfile[MAX_PATH+1] = "";
    char pdbfilepath[MAX_PATH+1] = "";
    Array2D <char> pdbids_2darray(MAX_NUM_PDBID, SIZE_PDBID+1);
    char ** pdbids = pdbids_2darray.array2D;
    bool isNonOptionArg = false;

    char *datadir = NULL;
    datadir = getenv("DATADIR");
    if (datadir == NULL || strcmp(datadir,"") == 0) {
        fprintf(stderr,"DATADIR not set, please specify the location for PDB files");
        return -1;
    } else {
        sprintf(pdbpath,"%s/%s",datadir,"pdb/data/structures/divided/pdb_dcp");
        sprintf(pdbobsoletepath,"%s/%s",datadir,"pdb/data/structures/obsolete/pdb_dcp");
        sprintf(pdbmodelspath,"%s/%s",datadir,"pdb/data/structures/models/current/pdb_dcp");
        sprintf(pdbmodelsobsoletepath,"%s/%s",datadir,"pdb/data/structures/models/obsolete/pdb_dcp");
    }

    i = 1;
    int cntPDBID = 0;
    while(i < argc )
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp(); 
                return 0;
            }
            else if( (strcmp(argv[i],"-l") == 0))  
            {
                my_strcpy(pdbidlistfile, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if( (strcmp(argv[i],"-d") == 0))  
            {
                my_strcpy(pdbpath, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if( (strcmp(argv[i],"-do") == 0))  
            {
                my_strcpy(pdbobsoletepath, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if( (strcmp(argv[i],"-o") == 0))  
            {
                my_strcpy(outfile, argv[i+1], MAX_PATH);
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
                fprintf(stderr,"Error! Wrong optional argument '%s'\n", argv[i]);
                return -1;
            }
        }
        else 
        {
            my_strcpy(pdbids[cntPDBID], argv[i], SIZE_PDBID);
            cntPDBID ++;
            if(cntPDBID > MAX_NUM_PDBID)
            {
                fprintf(stderr,"Error, number of pdbids > MAX_NUM_PDBID=%d\n", MAX_NUM_PDBID);
                return -1;
            }
            i ++;
        }
    }

    FILE *fpin = NULL;
    FILE *fpout = NULL;

    if( strcmp(pdbidlistfile,"") == 0 && cntPDBID == 0)
    {
        fpin=stdin; // read from stdin
    }
    else if (strcmp(pdbidlistfile, "") != 0)
    {
        fpin = fopen(pdbidlistfile,"r");
        checkfilestream(fpin,pdbidlistfile,"r"); 
    }
    if(fpin != NULL) cntPDBID = ReadInPDBIDs(fpin, pdbids, cntPDBID);

    if(fpin != stdin && fpin != NULL) fclose(fpin);

    if(strcmp(outfile, "") != 0) {
        fpout = fopen(outfile, "w");
        checkfilestream(fpout, outfile,"w");
    } else {
        fpout = stdout;
    }

    char *pstr;
    for(i = 0 ; i < cntPDBID; i++)
    {
        pstr = GetPDBFilePath(pdbids[i], pdbfilepath, pdbpath, pdbobsoletepath, pdbmodelspath, pdbmodelsobsoletepath);
        if(pstr != NULL) fprintf(fpout,"%s\n",pstr);
    }

    if(fpout != stdout && fpout != NULL) {
        fclose(fpout);
    }
    return EXIT_SUCCESS;
}
/*}}}*/
