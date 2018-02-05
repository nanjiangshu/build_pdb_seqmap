/*
 * =====================================================================================
 * 
 *        Filename:  splitfasta.cpp
 * 
 *     Description:  split a single fasta format file into individual fasta
 *                   files, each splited file consisting 1 sequence
 *           Usage:  splitfasta [options] fastafile
 * 
 *         Version:  1.0
 *         Created:  03/15/2006 06:00:06 PM CET
 *        Revision:  2006-11-02
 *        Compiler:  g++
 * 
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * 
 * =====================================================================================
 */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "myfunc.h"
#include "array.h"

#undef SIZE_ID
#define SIZE_ID MAX_PATH
int maxline = 10000;

void PrintHelp()
{
    fprintf(stdout,"Usage: splitfasta  [options] fastafile\n");
    fprintf(stdout,"Options:\n");
    fprintf(stdout,"\n");
    fprintf(stdout," -d <string>  : output to the defined directory, default = ./\n");
    fprintf(stdout," -e <string>  : extension for the output files, default = aa\n");
    fprintf(stdout," --maxline N  : set the maximual length of line in the fastafile, default = %d\n", maxline);
    fprintf(stdout," -h|--help    : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"split a dumped fasta file into a set of files including one sequence\n");
    fprintf(stdout,"each sequence is started by '>' as the head line\n");
    fprintf(stdout,"the name for the individual file is from the word following '>',\n");
    fprintf(stdout,"or rootname(fastafile)_cnt.aa if no word is following '>'\n");
    fprintf(stdout,"Create 2006-03-15, last modified 2007-05-22, Nanjiang Shu\n");
}
void SplitFasta(const char* seqDataFile, const char *outpath, const char *ext)/*{{{*/
{
	char seqFile[MAX_PATH+1] = "";
	//char command[MAX_COMMAND_LINE+1] = "";
	
	
	FILE *fpSeqData;
	FILE *fpSeq;
	fpSeqData = fopen(seqDataFile,"r");
    checkfilestream(fpSeqData, seqDataFile, "r");

	char id[SIZE_ID+1] = "";
	char str[SIZE_ID+10+1] = "";
    char rtname[MAX_PATH+1] = "";
    rootname(seqDataFile,rtname);
    int linesize;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 
    char delim[] = " \t\r\n>";
    char *pch = NULL;

	fpos_t pos;
    //int i;
	int cnt = 0 ;
	while((linesize = fgetline(fpSeqData,line,maxline)) != EOF)
	{
        if(linesize == 0 || IsBlankLine(line)) continue;

		if(line[0] == '>')//beginning of a record
		{
			cnt ++ ;
            my_strcpy(str, line, SIZE_ID+10);
            pch = strtok(str,delim);
            while(pch != NULL) //take the first word following '>' as id
            {
                my_strcpy(id, pch, SIZE_ID);
                break;
            }
            if(strcmp(id,"") == 0)// if no word following '>', set the id as rootname(seq-file)_cnt
            {
                sprintf(id,"%s_%d", rtname, cnt);
            }

			printf("%d\t:%s exported\n",cnt,id);
			sprintf(seqFile,"%s/%s.%s",outpath,id,ext);
			fpSeq = fopen(seqFile,"w");
            checkfilestream(fpSeq, seqFile, "r");
			fprintf(fpSeq,"%s\n",line);
			while(!feof(fpSeqData))
			{
				fgetpos(fpSeqData,&pos);
				linesize = fgetline(fpSeqData,line,maxline);
                if(linesize == 0 || IsBlankLine(line)) continue;
				if(line[0] == '>')
				{
					fsetpos(fpSeqData,&pos);
					break;
				}
				else
					fprintf(fpSeq,"%s\n",line);

			}
			fclose(fpSeq);
		}
	}
	fclose(fpSeqData);
	printf("%d sequences exported!\n",cnt);
}
/*}}}*/
int main(int argc, char** argv)/*{{{*/
{
    int i;
    char outpath[MAX_PATH+1] = "";
    char fastafile[MAX_PATH+1] = "";
    char ext[MAX_PATH+1] = "aa";

	char *returnGetCWD=NULL;
    
    returnGetCWD = getcwd(outpath,MAX_PATH); // get the current outpath

    if(argc < 2) 
    {
        PrintHelp();
        return -1;
    }

    i = 1;

    while(i < argc)
    {
        if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
        {
            PrintHelp(); 
            return 0;
        }
        else if( (strcmp(argv[i],"-d") == 0))  
        {
            my_strcpy(outpath,argv[i+1],MAX_PATH);
            VerifyFolder(outpath);
            i += 2;
        }
        else if( (strcmp(argv[i],"-e") == 0))  
        {
            my_strcpy(ext,argv[i+1],MAX_PATH);
            i += 2;
        }
        else if( (strcmp(argv[i],"--maxline") == 0))  
        {
            maxline = atoi(argv[i+1]);
            i += 2;
        }
        else 
        {
            my_strcpy(fastafile,argv[i], MAX_PATH);
            i += 1;
        }
    }

    if(strcmp(fastafile,"") == 0)
    {
        fprintf(stderr,"Error! fastafile not set\n");
        return -1;
    }

    SplitFasta(fastafile, outpath,ext);

    return (0);
}
/*}}}*/
