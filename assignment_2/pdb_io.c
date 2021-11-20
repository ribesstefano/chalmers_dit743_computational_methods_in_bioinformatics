#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_LENGTH	81
	/*
	 * See http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/part_10.html
	 * There are 80 characters on a line.  Plus the newline character = 81.
	 */


int main(argc, argv)
	int	argc;
	char	**argv;
{
	FILE	*stream;
	char	line[LINE_LENGTH];


	/*
	 * Declare strings to hold parts of an ATOM record.
	 * See http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/part_62.html
	 *
	 * Each string is declared to have a length that is one greater than
	 * the field size.  This is to accommodate a null character.
	 */

	char	s_serial[6];
	char	s_name[5];	/* Alpha-carbon is " CA "; calcium is "CA  " */
	char	s_altLoc[2];	/* Usually " " */
	char	s_resName[4];
	char	s_chainID[2];
	char	s_resSeq[5];
	char	s_iCode[2];	/* Usually " " */
	char	s_x[9];
	char	s_y[9];
	char	s_z[9];

	int	serial;
	int	resSeq;
	double	x;
	double	y;
	double	z;

	if ( argc<2 ) {
		(void) fprintf(stderr, "usage: pdb_io file\n");
		exit(0);
	}

	if ( (stream = fopen(argv[1], "r")) == NULL ) {
		(void) fprintf(stderr, "Unable to open %s\n", argv[1]);
		exit(0);
	}

	while ( fgets(line, LINE_LENGTH, stream) ) {
		if ( strncmp(line, "ATOM  ", 6) == 0 ) {

			/*
			 * Split the line into its constituent fields.
			 * We are only interested in columns 1-54.
			 */

			strncpy(s_serial,  &line[6],  5); s_serial[5]  = '\0';
			strncpy(s_name,    &line[12], 4); s_name[4]    = '\0';
			strncpy(s_altLoc,  &line[16], 1); s_altLoc[1]  = '\0';
			strncpy(s_resName, &line[17], 3); s_resName[3] = '\0';
			strncpy(s_chainID, &line[21], 1); s_chainID[1] = '\0';
			strncpy(s_resSeq,  &line[22], 4); s_resSeq[4]  = '\0';
			strncpy(s_iCode,   &line[26], 1); s_iCode[1]   = '\0';
			strncpy(s_x,       &line[30], 8); s_x[8]       = '\0';
			strncpy(s_y,       &line[38], 8); s_y[8]       = '\0';
			strncpy(s_z,       &line[46], 8); s_z[8]       = '\0';

			/*
			 * Print the string values of the fields.
			 */

			// printf("serial= \"%s\"\n", s_serial);
			// printf("name=   \"%s\"\n", s_name);
			// printf("altLoc= \"%s\"\n", s_altLoc);
			// printf("resName=\"%s\"\n", s_resName);
			// printf("chainID=\"%s\"\n", s_chainID);
			// printf("resSeq= \"%s\"\n", s_resSeq);
			// printf("iCode=  \"%s\"\n", s_iCode);
			// printf("x=      \"%s\"\n", s_x);
			// printf("y=      \"%s\"\n", s_y);
			// printf("z=      \"%s\"\n", s_z);

			/*
			 * Convert the numeric fields to integers or doubles.
			 * The library functions atoi() and atof() are
			 * described in the UNIX manual pages ('man atoi' and
			 * 'man atof').
			 */

			serial = atoi(s_serial);
			resSeq = atoi(s_resSeq);
			x      = atof(s_x);
			y      = atof(s_y);
			z      = atof(s_z);

			/*
			 * Print the line that was read, then use printf() to
			 * print the individual fields in Protein Data Bank
			 * format.
			 * We are only interested in columns 1-54.
			 * The format of printf's format string is described in
			 * the UNIX manual pages ('man 3 printf').
			 */

			printf("%s\n", line);

			printf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f\n",
				serial,
				s_name,
				s_altLoc,
				s_resName,
				s_chainID,
				resSeq,
				s_iCode,
				x,
				y,
				z);
		}
	}


        return 0;
}
