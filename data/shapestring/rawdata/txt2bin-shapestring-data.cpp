/*
 * =====================================================================================
 *
 *       Filename:  txt2bin-shapestring-data.cpp
 *
 *    Description:  change the shape string raw data, shape string region
 *                  definition on the Ramachandran plot from ascii format to
 *                  binary format
 *
 *        Version:  1.0
 *        Created:  11/05/2007 04:15:51 PM CET
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 *
 * =====================================================================================
 */

#include <cstdio>
#include <cstring>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"

/*ChangeLog 2010-04-14 
 * The datatype for PointDef changed to unit8 from int
 * */

struct PointDef
{
     unit8 x;
     unit8 y;
     unit8 intensity;
};

void PrintHelp()
{
    fprintf(stdout,"Usage: txt2bin-shapestring-data [options] rawdata-file\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  -o|--out <outfile>  : output the result to outfile, default = rootname.bin\n");
    fprintf(stdout,"  -h|--help           : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-11-05, updated 2010-04-14, Nanjiang Shu\n");
}
void PrintVerboseHelp() { }

int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;
    if(argc < 2) 
    {
        PrintHelp();
        return 0;
    }
    int i,j;
    char infile[MAX_PATH+1] = "";
    char outfile[MAX_PATH+1] = "";
    double value = 0.0;
    const char control_option[] = ""; //options which control the program, and does not take parameters
    bool isAll = false;
    bool isQuiet = false;
    bool isSingle = false;

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                i ++;
            }
            else if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp(); 
                return 0;
            }
            else if(strcmp(argv[i],"-H") == 0 )
            {
                PrintVerboseHelp();
                return 0;
            }
            else if( (strcmp(argv[i],"-o") == 0) || (strcmp(argv[i], "--out") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, outfile)) == -1)
                    return -1;
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
        else //non-option argument
        {
            my_strcpy(infile, argv[i], MAX_PATH);
            i ++;
        }
    }/*}}}*/


    FILE *fpout = NULL;
    if(strcmp(outfile,"") == 0 || strcasecmp(outfile, "stdout") == 0)
    {
        fpout = stdout;
    }
    else
    {
        fpout = fopen(outfile, "w");
        checkfilestream(fpout, outfile,"w");
    }

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    Array1D <PointDef> points_1darray(32401);
    PointDef *points = points_1darray.array1D;
    FILE *fpin = fopen(infile, "r");
    if(checkfilestream(fpin, infile, "r") == -1) return -1;
    int sscanf_status = 0;
    i = 0;
    int a,b,c;
    while((linesize = fgetline(fpin, line, maxline )) != EOF)
    {
        sscanf_status = sscanf(line,"%d %d %d", &a, &b, &c);
        assert(sscanf_status == 3);
        points[i].x = unit8(a); 
        points[i].y = unit8(b); 
        points[i].intensity = unit8(c);
        i ++;
    }
    assert (i <= 32400);
    fclose(fpin);
    int num = i;

    //write to binary file
    char rtname[MAX_PATH+1] = "";
    rootname(infile, rtname);
    char binaryfile[MAX_PATH+1] = "";
    sprintf(binaryfile,"%s.bin", rtname);
    FILE *fpwb = fopen(binaryfile, "wb");
    if(checkfilestream(fpwb, binaryfile, "wb") == -1) return -1;
    fwrite(&num, sizeof(int), 1, fpwb);
    fwrite(points, sizeof(PointDef), num, fpwb);
    fclose(fpwb);
    fprintf(stdout,"%s -> %s\n", infile, binaryfile);;


    if(fpout != NULL && fpout != stdout) fclose(fpout);

    return 0;
}
/*}}}*/
