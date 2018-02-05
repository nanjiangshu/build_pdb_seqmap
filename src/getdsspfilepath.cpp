#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include "myfunc.h"

#ifndef MAX_PATH
#define MAX_PATH 266
#endif

#undef NN
#define NN 100
/*****************************************************************************
 * ChangeLog 2008-02-06
 *  add the option: -l pdbidListFile
 ****************************************************************************/

void PrintHelp()
{
    fprintf(stdout, "Usage: getdsspfilepath [options] pdbid1 pdbid2 ...\n");
    fprintf(stdout, "    pdbid is case insensitive\n");
    fprintf(stdout, "options:\n");
    fprintf(stdout, "   -d path          : path for dssp files, default = $DATADIR/dssp\n");
    fprintf(stdout, "   -l pdbidListFile : supply the pdbid list file\n");
    fprintf(stdout, "   -h|--help        : print this help message\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Created 2006-01-01, updated 2010-01-18, Nanjiang Shu\n");
}

int main(int argc , char** argv)/*{{{*/
{

    if ( argc < 2)
    {
        PrintHelp();
        return 1;
    }

	int i;
    char dssppath[MAX_PATH+1] = "";
	char dsspfilepath[MAX_PATH+1] = "";
    char idListFile[MAX_PATH+1] = "";
    char *pstr;
    char pdbid[SIZE_PDBID+1] = "";

    char *datadir             = NULL;
    datadir = getenv("DATADIR");
    if (datadir == NULL || strcmp(datadir,"") == 0)
    {
        strcpy(dssppath, "/data/dssp");
    }
    else
    {
        sprintf(dssppath,"%s/%s",datadir,"dssp");
    }

    int cnt_id_arg = 0; /*count the number of pdbids supplied in the argument list*/
    bool isNonOptionArg = false;
    i = 1;
	while(i < argc )/*{{{*/
	{
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 ) { PrintHelp(); break; }
            else if(strcmp(argv[i],"-d") == 0)
            {
                my_strcpy(dssppath,argv[i+1],MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"-l") == 0)
            {
                my_strcpy(idListFile,argv[i+1], MAX_PATH);
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
            my_strcpy(pdbid, argv[i], SIZE_PDBID);
            my_strlwr(pdbid);
            pstr = GetDSSPFilePath(pdbid,dsspfilepath,dssppath);
            if(pstr != NULL){ fprintf(stdout, "%s\n",pstr);}
            i += 1;
            cnt_id_arg ++;
        }
	}/*}}}*/

    char id[NN+1] = "";
    int strsize= 0;
    if(cnt_id_arg == 0)// if no id has been set in the argument list, read from stdin
    {
        FILE *fpin = NULL;
        if(strcmp(idListFile,"") != 0)
        {
            fpin = fopen(idListFile, "r");
        }
        else
        {
            fpin = stdin;
        }

        if(checkfilestream(fpin, idListFile, "r") != -1)
        {
            while((strsize=fgetdelim(fpin, id, WHITE_SPACE, NN)) != EOF)
            {
                if (strsize < 4)
                {
                    continue;
                }
                my_strcpy(pdbid, id, SIZE_PDBID);
                my_strlwr(pdbid);
                pstr = GetDSSPFilePath(pdbid,dsspfilepath,dssppath);
                if(pstr != NULL) { fprintf(stdout, "%s\n",pstr); }
            }
        }
        if(fpin != stdin && fpin != NULL)
        {
            fclose(fpin);
        }
    }
	return 0;
}
/*}}}*/

