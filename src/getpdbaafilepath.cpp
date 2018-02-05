#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "mypro.h"
#include "myfunc.h"

//ChangeLog/*{{{*/
/*****************************************************************************
 * ChangeLog 2007-10-10
 *  add --idtype, so that when idtype = 1, pdbaafilename = $pdbaapath/id.aa
 * ChangeLog 2008-02-06
 *  add the option: -l idListFile
 ****************************************************************************/
/*}}}*/

void PrintHelp()
{
    fprintf(stdout,"Usage: getpdbaafilepath [options] id1 id2 ...\n");
    fprintf(stdout,"  id is the standardized pdb chain identifier, case insensitive\n");
    fprintf(stdout,"  the first 6 characters of an id is used if there are more\n");
    fprintf(stdout,"  piping supported\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  -d path      : path for pdbaa files, default = $DATADIR/pdbaa\n");
    fprintf(stdout,"  --idtype 0|1 : 0 -- standardized id, 1 -- exact id name, default=0\n");
    fprintf(stdout,"  -l idListFile: supply the idList file\n");
    fprintf(stdout,"  -h|--help    : print this help message\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created 2006-02-04, updated 2010-04-15, Nanjiang Shu\n");
}

#define NN 30
int main(int argc , char** argv)/*{{{*/
{
    int i;
    char pdbaapath[MAX_PATH+1] = "";
    char idListFile[MAX_PATH+1] = "";
    int idtype = 0;
    char pdbaafilepath[MAX_PATH+1] = "";
    char *pstr;
    char id[NN+1] = "";


    char *datadir = NULL;
    datadir = getenv("DATADIR");
    if (datadir == NULL || strcmp(datadir,"") == 0)
    {
        strcpy(pdbaapath, "/data/pdbaa");
    }
    else
    {
        sprintf(pdbaapath,"%s/%s",datadir,"pdbaa");
    }

    int cnt_id_arg =0;
    bool isNonOptionArg = false;

    i = 1;

    while(i < argc )/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp();
                return 0;
            }
            else if(strcmp(argv[i],"-d") == 0)
            {
                my_strcpy(pdbaapath,argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"-l") == 0)
            {
                my_strcpy(idListFile,argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--idtype") == 0)
            {
                idtype=atoi(argv[i+1]);
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
            my_strcpy(id, argv[i], NN);
            if(idtype == 0)
            {
                StdID(id);
                //LOG: 2007-04-25 12:19:15 Wednesday  Week 16 <nanjiang@casio.fos.su.se>
                // first convert the id to StdID, since GetPDBAAFilePath accept standard chain
                // id
                pstr = GetPDBAAFilePath(id,pdbaafilepath,pdbaapath);
                if(pstr != NULL)
                {
                    fprintf(stdout,"%s\n",pstr);
                }
            }
            else { fprintf(stdout,"%s/%s.aa\n", pdbaapath, id); }
            cnt_id_arg ++;
            i += 1;
        }
    }/*}}}*/

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
            while( fgetdelim(fpin, id, WHITE_SPACE, NN) != EOF)
            {
                StdID(id);
                if(idtype == 0)
                {
                    pstr = GetPDBAAFilePath(id,pdbaafilepath,pdbaapath);
                    if(pstr != NULL)
                        fprintf(stdout,"%s\n",pstr);
                }
                else
                {
                    fprintf(stdout,"%s/%s.aa\n", pdbaapath, id);
                }
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

