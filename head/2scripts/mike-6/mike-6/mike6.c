#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>

void load_multi_fasta(char **Seq,FILE *fic, char **name, int *lg , int *nb_seq);
char ** taballoc(int nl,int nc);
void tabdesalloc ( char **matrice,int nl);
char * load_single_fasta(FILE *fic, char *name, int *lg);
double ** ftaballoc(int nl,int nc);
int ** itaballoc(int nl,int nc);
int is_nt (char nt);
char genetic_code_hmtDNA (char nt1, char nt2, char nt3);
char complement (char c);
double ** order (double *tab, int nb);
double ** rank( double **mat1, int **mat2, int nb);
int ** order_int (int *tab, int nb);
double spearman (double ** rang, int nb);
double Pearson (double *LD, int *dist, int nb);
double permute_spearman (double ** rang, int permut, double rho, int nb);
double permute_Pearson (double *LD, int *pos, int permut, double rho, int nb);
double permute_spearman2 (double ** LDordonne, int *pos, int permut, double rho, int nb) ;
static void PrintDate ();



/* Constantes definies pour ran 1 */
#define IA 16807     
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
//double ran1 (long *idum); fonction remplacee par Rand

#define MAXLENGTH 10000
#define MAXSEQ 190

void 				SetSeed (int seed);
double 				rndu (void);
static int				 z_rndu=137;
static unsigned 		w_rndu=13757;
static int 			Rand(int max);


int main()
{
FILE *fich_seq, *fich_name,*fich_seq_ref, *fich_pos, *fich_pos_syn, *fich_doublets;  /* fichier sequence en entree format FASTA */
FILE *fich_outnames, *fich_out_syn;
char **seq, **haplo,**name, nameref[100],c,aa[3],bb, Min1,Maj1,Min2, Maj2,nt1,nt2,nt3;
char *nom_fich_seq, **nom_fich, *nom_fich_pos, *nom_fich_syn,*all_maj,*all_min, *name1, *nom_doublets,*buffer;
int compteurH,nbhaplo,Laln,nb_fich,lg_seq[100], lg_seqref, nbligne, nb_lg, nb_cl,i, j, kk, z, nbseq, nbpos, nb_site_pol, nb_site_pol_syn,length, nb_doublets;
int zzzz,*typepol,*haplo_eff,*pos,*pos_syn, *pos_doublets, *distance, *test2all, **mat2,zdoublets;
int  compteur,compteur2,compteur3, cMin1Min2, cMin1Maj2, cMaj1Min2, cMaj1Maj2,min, *pos_nr, *posi_nr, nb_pos_syn, zz, nbpermut,ts,tv;
double testDprime,*Min, *Dprime, *D, *r_square, Dmax, **mat1, **mat1b,**rang,**rangb,rho_pearson_r2, rho_pearson_D,rho_spearman_r2,rho_spearman_D,P;
//long idum;
long loc;

char 	*infilename;

infilename = (char*) calloc (200, sizeof(char));

/*ouverture des fichiers de sequence*/
printf("enter your file's name>>");
scanf("%s",infilename);

printf("\nRelationship between LD and distance - Mike6.c Version 1.00\n");
printf("______________________________________________________________________________________________\n\n");
PrintDate ();

printf("\nInput file: %s\n",infilename);

nom_fich_syn=calloc(60,sizeof (char));
loc=0;
SetSeed(time(NULL));

//printf("enter a negative integer to initialize ran1 :");
//scanf("%d",&idum); suppression de l'entree d'un nb aleatoire negatif par l'utilisateur, car remplacement  de la fct ran1 par Rand

/* allocation memoire */

if( (seq = taballoc (MAXSEQ,MAXLENGTH)) == NULL)
	{
        fprintf(stderr,"\nERROR: not enough memory (calloc).\n");
        exit(1);
       	}
if( (haplo = taballoc (MAXSEQ,MAXLENGTH)) == NULL)
	{
        fprintf(stderr,"\nERROR: not enough memory (calloc).\n");
        exit(1);
       	 }
buffer=calloc(MAXLENGTH,sizeof (char));

fich_seq = fopen(infilename,"r");

if (fich_seq == NULL) {printf("fich_seq unopened\n"); exit(1);}
load_multi_fasta(seq,fich_seq, name, lg_seq, &nbseq);

fclose(fich_seq);
for (i=0;i<(nbseq-1);i++)
	{
	if (lg_seq[0]!= lg_seq[i+1]) { printf("pb: not the same number of nucleotides in sequences %d : %d first seg %d",i+1,lg_seq[i+1],lg_seq[0]); exit(1);}
	}
nbpos=lg_seq[0];

test2all=calloc( nbpos, sizeof(int));

/*verifie qu'il n'y a que 2 alleles segregeant sur un site*/

nb_site_pol=0;
all_maj=calloc(nbpos,sizeof(char));
all_min=calloc(nbpos,sizeof(char));
Min=calloc(nbpos, sizeof(double));	
pos_nr=calloc(nbpos, sizeof(int));
posi_nr=calloc(nbpos, sizeof(int));
typepol=calloc(nbpos, sizeof(int));
haplo_eff=calloc(nbseq, sizeof(int));
zz=0;
zdoublets=0;
min=0;
ts=tv=0;

/*1ere lecture
printf("\n protein analysed... \n");*/

Laln=nbpos;
for (j=0;j<nbpos;j++)
	{
	compteur=0;
	
	for (i=0;i<nbseq;i++)
		{
		
		if (seq[i][j]!='a' && seq[i][j]!='t' && seq[i][j]!='g' && seq[i][j]!='c' )
			{
			++compteur;
			}
		}
	if (compteur!=0) --Laln;
	}
//printf(" \nLaln %d", Laln); //nb de nt alignes
//printf("\n<b>More than 2 alleles segregating at position :</b>\n");
compteur3=0;

for (j=0;j<nbpos;j++)
	{
	compteur=0;
	compteur2=1;
	test2all[j]=0;
	zz=0;
	for (i=1;i<nbseq;i++)
		{

			if ( seq[0][j] != seq[i][j])  
				{ 
					if (compteur!=0)
						{
							aa[2]=seq[i][j];

							if (aa[1] != aa[2])
								{
								if (zz==0) 
									{
									compteur3++;
									if(compteur3==1)
										{
										 printf("\nMore than 2 alleles segregating at position :\n");
										 }

									//printf("\nmore than 2 nt segregating pos %d\t%c\t%c\t%c",j,aa[1],aa[2],seq[0][j]);}
									printf("%d - ", j);
									}
									zz=1;
									test2all[j]=1;
								}
						}
					else
						{
							aa[1]=seq[i][j];
							++compteur;
						}
					}
					else {++compteur2;}
		}


			if ( test2all[j]==0 &&compteur2!=nbseq /* compteur2<(nbseq-1) */  && is_nt(seq[0][j])==1 && is_nt(aa[1])==1) /*restriction to informative sites*/

				{
					min= (compteur2 < (nbseq-compteur2) ) ? compteur2 : (nbseq-compteur2) ;
					all_min[nb_site_pol]= (compteur2 < (nbseq-compteur2) ) ? seq[0][j] : aa[1] ;
					all_maj[nb_site_pol]= (compteur2 < (nbseq-compteur2) ) ? aa[1] : seq[0][j] ;
					posi_nr[nb_site_pol]=j;
					Min[nb_site_pol]= (float) (min)/ (float) (nbseq);


					if (j%3==0) //nt1 pol
						{
						nt1=aa[1];
						nt2=seq[0][j+1];
						nt3=seq[0][j+2];
						
						if (genetic_code_hmtDNA(nt1,nt2,nt3)==genetic_code_hmtDNA(seq[0][j],nt2,nt3)) 
								{
								typepol[nb_site_pol]=1;

								
								}
						else		
								{
								typepol[nb_site_pol]=0;

								
								}
						
						}
					else if ((j+1)%3==0) //nt3
						{
						nt1=seq[0][j-2];
						nt2=seq[0][j-1];
						nt3=aa[1];
						
						if (genetic_code_hmtDNA(nt1,nt2,nt3)==genetic_code_hmtDNA(nt1,nt2,seq[0][j])) 
								{
								typepol[nb_site_pol]=1;

								}
						else		
								{
								typepol[nb_site_pol]=0;

								}
						
						}
					else //nt2
						{
						nt1=seq[0][j-1];
						nt2=aa[1];
						nt3=seq[0][j+1];
						if (genetic_code_hmtDNA(nt1,nt2,nt3)==genetic_code_hmtDNA(nt1,seq[0][j],nt3)) 
								{
								typepol[nb_site_pol]=1;

								}
						else		
								{
								typepol[nb_site_pol]=0;
								}
						}
						
					++nb_site_pol;

				}
	}
if(compteur3>0) {printf("\nThese positions were removed from LD-distance analysis\n");}

printf("\n\nNumber of sequences:%d \nNumber of aligned nucleotides :%d \n\nNumber of polymorphic sites:%d\n",nbseq,nbpos,nb_site_pol);

if(nb_site_pol<4){
				printf("\n------------------------------  TOO FEW SITES TO ESTIMATE LD - DISTANCE CORRELATION COEFFICIENT------------------------------ \n");
				exit(1);
			      }
//get haplotypes
for (j=0;j<lg_seq[0];j++)
	{
	haplo[0][j]=seq[0][j];
	}
nbhaplo=1;
haplo_eff[0]=1;
for (i=1;i<nbseq;i++)
	{

	haplo_eff[i]=0;
	compteurH=0;
	for (z=1;z<=nbhaplo;z++)
		{
		compteur=0;
		for (j=0;j<lg_seq[0];j++)
			{
			if (seq[i][j]!=haplo[z-1][j]) //&& seq[i][j]!='-'
				{
					++compteur; //if (is_nt(seq[i][j])==0) {printf("ind %d pos %d %c\n",i,j,seq[i][j]);}
				}
			}
		if (compteur==0) //same haplotype
			{
					++haplo_eff[z-1];
					//printf("same haplotype haplo_eff[%d]=%d\n",z-1,haplo_eff[z-1]);
			}
		else
	       		{
					compteurH++;

				}
		}
	if (compteurH==nbhaplo) //there is a new haplotype
		{

		for (j=0;j<lg_seq[0];j++) 
			{
				haplo[z-1][j]=seq[i][j];
			}
		haplo_eff[z-1]=1;
		nbhaplo++;
		}
	}

/*distance calculation*/

compteur=nb_site_pol*(nb_site_pol-1)/2;

distance= calloc( compteur+1, sizeof(int));
if (distance== NULL)
	 {printf("\ndistance not allocated\n");exit(1); }

Dprime=calloc(compteur+1,sizeof(double));
D=calloc(compteur+1,sizeof(double));
r_square=calloc(compteur+1,sizeof(double));
if (D == NULL) {printf("\nD not allocated\n"); exit(1);}
if (Dprime == NULL) {printf("\nDprime not allocated\n");exit(1);}
if (r_square == NULL) {printf("\nr_square not allocated\n");exit(1);}

compteur=0;
testDprime=0;
for (i=0;i<(nb_site_pol-1);i++)
	{

	Min1=all_min[i];
	Maj1=all_maj[i];
	for (j=(i+1);j<nb_site_pol;j++)
		{

			Min2=all_min[j];
			Maj2=all_maj[j];
			distance[compteur]= ( abs(posi_nr[i]-posi_nr[j]) < (16569-abs(posi_nr[i]-posi_nr[j])) ) ? abs(posi_nr[i]-posi_nr[j]) : (16569-abs(posi_nr[i]-posi_nr[j]));

			cMin1Min2=0;
			cMin1Maj2=0;
			cMaj1Min2=0;
			cMaj1Maj2=0;		
			for (kk=0;kk<nbseq;kk++)
				{
					if (seq[kk][posi_nr[j]] == Min2 && seq[kk][posi_nr[i]] == Min1) 
						{
							++cMin1Min2;

						}
					if (seq[kk][posi_nr[j]] == Min2 && seq[kk][posi_nr[i]] == Maj1) 
						{
							++cMaj1Min2;

						}	
					if (seq[kk][posi_nr[j]] == Maj2 && seq[kk][posi_nr[i]] == Maj1) 
						{
							++cMaj1Maj2;

						}	
					if (seq[kk][posi_nr[j]] == Maj2 && seq[kk][posi_nr[i]] == Min1) 
						{
							++cMin1Maj2;

						}							
				}
		
			D[compteur]= (float) (cMaj1Maj2*cMin1Min2)/ (float) (nbseq*nbseq) - (float) (cMin1Maj2*cMaj1Min2)/ (float) (nbseq*nbseq);

			Dmax=1;
			if (D[compteur]>0)
				{
					Dmax= ( (Min[i] * (1-Min[j])) < (Min[j]* (1-Min[i])) ) ? (Min[i] * (1-Min[j])) : (Min[j]* (1-Min[i]));

				}
			
		
			if (D[compteur]<0) 
				{
					Dmax= ( (Min[i] * Min[j]) < ((1-Min[j])*(1-Min[i])) ) ? (Min[i] * Min[j]) : ((1-Min[j])*(1-Min[i]));

				}	
				
			Dprime[compteur]=fabs( D[compteur]/Dmax );
			testDprime+=Dprime[compteur];
			r_square[compteur]=  D[compteur]*D[compteur]/ (Min[i] * (1-Min[i])* Min[j]*(1-Min[j])) ;
			compteur++;

		}
	}
mat1=ftaballoc(compteur,2);
if (mat1 == NULL) 
	{
	printf("\nnot enough memory to allocate mat1\n");
	exit(1);
	}

mat1b=ftaballoc(compteur,2);
if (mat1b == NULL)
	{
	printf("\nnot enough memory to allocate mat1b\n");
	exit(1);
	}
	
mat2=itaballoc(compteur,2);
if (mat2 == NULL) 
	{
	printf("\nnot enough memory to allocate mat2\n");
	exit(1);
	}
	
rang=ftaballoc(compteur,2);
if (rang == NULL) 
	{
	printf("\nnot enough memory to allocate rang\n");
	exit(1);
	}
	
rangb=ftaballoc(compteur,2);
if (rangb == NULL) 
	{
	printf("\nnot enough memory to allocate rangb\n");
	exit(1);
	}
mat1 = order (r_square, compteur);
mat2 = order_int (distance,compteur);
/*for (i=0; i<compteur; i++)
	{
	printf ("%i %f \t %f %f\t\t %f\n",i, r_square[i], mat1[i][1], mat1[i][0],Dprime[i]);

	}*/

	
rang= rank(mat1,mat2, compteur);
if (rang == NULL) {printf("\npb ds la fonction rank : le tableau de sortie est vide!!\n");
					exit(1);}
rho_pearson_r2=Pearson(r_square,distance,compteur);
printf("\n\n---- PEARSON CORRELATION COEFFICIENT BETWEEN LINKAGE DESEQUILIBRIUM MEASURED AS R SQUARE AND DISTANCE ----\n");
//printf("\nPearson correlation coefficient between LD measured as r square and distance = %f\n", rho_pearson_r2);

if (nb_site_pol<193) nbpermut=1000;
else  nbpermut=500;
P=permute_Pearson (r_square, posi_nr, nbpermut, rho_pearson_r2, nb_site_pol);
//printf("\rP-value(LD r square_d) =  %f\t (number of randomizations : %d)\n",P,nbpermut);

printf("\n\tCorrelation coefficient : %.3f\n\tP-value : %.3f\n\tNumber of randomizations : %d\n\n" ,rho_pearson_r2,P,nbpermut);

if (testDprime<(compteur-0.1))
	{
		rho_pearson_D=Pearson(Dprime,distance,compteur);
		printf("\n---- PEARSON CORRELATION COEFFICIENT BETWEEN LINKAGE DESEQUILIBRIUM MEASURED AS |D'| AND DISTANCE----\n");
		//printf("\n\nPearson correlation coefficient between LD measured as |D'| and distance = %f\n", rho_pearson_D);

		mat1b = order (Dprime, compteur);	
		rangb= rank(mat1b,mat2, compteur);
		if (rangb == NULL) {printf("\npb ds la fonction rank : le tableau de sortie est vide!!\n");}

	}
else
	{
		//printf("\nDprime is always 1\t");
		printf("\n|D'| ------------------------------ IS EQUAL TO ONE FOR MOST PAIRS OF SITES------------------------------ \n------------------------------ NOT ENOUGH DATA TO ESTIMATE PEARSON CORRELATION COEFFICIENT WITH LD MEASURED AS |D'|------------------------------ ");
		//printf("D'1\t");
	}
if (testDprime<(compteur-0.1))
	{
		if (nb_site_pol<193) nbpermut=1000;
		else  nbpermut=500;
		P=permute_Pearson (Dprime, posi_nr, nbpermut, rho_pearson_D, nb_site_pol);
		//printf("\rP-value(LD |D'|_d) =  %f\t (number of randomizations : %d)\n",P, nbpermut);

		printf("\n\tCorrelation coefficient : %.3f\n\tP-value : %.3f\n\tNumber of randomizations : %d\n\n" ,rho_pearson_D,P,nbpermut);
	}

free(distance);
free(Min);

free(test2all);
free(all_maj);
free (all_min);
free(pos_nr);
free(posi_nr);
free(typepol);
free(D);
free(Dprime);
free(r_square);
for (j=0;j<compteur;j++)
	{
	free(mat1b[j]);
	free(mat2[j]);
	}


return 0;
}

/*****************  load_multi_fasta ***************************/
void load_multi_fasta(char **SS, FILE *fic, char **name, int *lg, int *nb_seq)
{
char string[100], name1[50];
int i,j=0, length = 0, maxlength=0,z;
long pos;
char c;
/* 1ere passe */
fgets(string, 100, fic);
if(string[0] != '>') {printf("\nformat error ! not a fasta file %c\n",string[0]);}
if( sscanf(string+1, "%s", name1) == EOF) {
        printf("\nEXIT: Format error ! Not a FASTA file\n %d %s\n", sscanf(string+1, "%s", name1), string+1);
        exit(1);
       		 }
 
pos = ftell(fic);	
  rewind(fic); 
  fscanf(fic,"%s\n",name1);
i=0;

while((c = fgetc(fic)) != EOF) 
	{
		c=tolower(c);
		
		/**/
		
		if (c!= '>')
		{
        if(c != 'a' && c != 't' && c != 'g' && c!='c' && c!='-' && c!='y' && c!='w' && c!='s' && c != '?' && c!='n' && c!='k' && c!='r' && c!='m' && c!='b'  && c!='h')
         	{
         	if(c == ' ' || c == '\n' || c == '\t' || c == '\r' ) continue;
        	else 
        		{printf("\n[%c] not allowed caracter\n",c);
        		exit(1);
        		}
         	}
        
        else  {length++;
        		/*printf("%c",c);*/}
        }
        
        else
        {
        lg[j]=length;
        
        ++j;
        
        maxlength = (length > maxlength ? length : maxlength);
        
        length=0;

        fscanf(fic,"%s\n",name1);
	
        pos = ftell(fic);
        
		
        }
        	
    }
 
 lg[j]=length;

        ++j;
        
        maxlength = (length > maxlength ? length : maxlength);
        
        length=0;
 
* nb_seq = j;
      
if (j>MAXSEQ) {printf("%d is greater than maximum authorisez nb of seg %d\n",j,MAXSEQ); exit(1);}
if(maxlength>MAXLENGTH) { printf("%d is > than maximum authorised length %d\n",maxlength,MAXLENGTH); exit(1);}
/*if( (name = taballoc (j+1,50)) == NULL) {
        fprintf(stderr,"\nERROR: not enough memory (calloc).\n");
        exit(1);
        }*/
for (i=0;i<j;i++)
	{
	for (z=0;z<maxlength;z++)
		{SS[i][z]='0';
		}
	}
/* charge la sequence */
pos=0;
fseek(fic, pos, 0);
z = 0;
i=0;
c = fgetc(fic);
if ( c != '>')
			{
       		printf("error first char not '>' !!!\n");
  			}
fscanf(fic,"%s\n",name1);
/*name[i] = string;*/
pos = ftell(fic);
  			 			
while((c = fgetc(fic)) != EOF ) 
	{
		c=tolower(c);
		
        if (c == 'a' || c =='t' || c =='g' || c =='c' || c =='-' || c=='y' || c=='w' || c=='s' || c=='n' || c=='k' || c=='r' || c=='m' || c=='b')
        												{
        												if (z>=maxlength)
        														{
        														printf("\n error ! more tham maxlength caracter %i\n",i);
        														exit(1);
        														}
        												SS[i][z] = c;
        											
        												z++;
        												
        											
        												} 
      		     
        else if (c == '>')
        	 {
        	if (z !=lg[i] || z!=maxlength) {printf("\n ? seq %d lg %d 2nd lg %d maxlength %d\n",i,lg[i],z,maxlength);}
        	z=0;
        	i++;
        	fscanf(fic,"%s\n",name1);
        	/*name[i]=string;*/
			pos = ftell(fic);
			
        	}
        
     }
	
/*printf("%d sequences of maxlength %d entered\n", i+1, maxlength);*/
}
/***************** end load_multi_fasta ***************************/
char ** taballoc(int nl,int nc)
{
char **pmat;
int q;
pmat=(char **)calloc(nl,sizeof(char *)); 
for (q=0;q<nl;q++) *(pmat+q)=(char *)calloc(nc,sizeof(char));
return pmat;
}
int ** itaballoc(int nl,int nc)
{
int **pmat;
int q;
pmat=(int **)calloc(nl,sizeof(int *)); 
for (q=0;q<nl;q++) *(pmat+q)=(int *)calloc(nc,sizeof(int));
return pmat;
}
double ** ftaballoc(int nl,int nc)
{
double **pmat;
int q;
pmat=(double **)calloc(nl,sizeof(double *)); 
for (q=0;q<nl;q++) *(pmat+q)=(double *)calloc(nc,sizeof(double));
return pmat;
}
void tabdesalloc( char **matrice,int nl)
{
int q;
for (q=0;q<nl;q++) free(matrice+q);
printf("espace memoire desalloue \n");
}
char complement (char c)
	{
	char nt;
	
	if (c== 'a') nt='t';
	if (c=='t') nt='a';
	if (c=='g') nt='c';
	if (c=='c') nt='g';
	
	if (nt != 'a' && nt!='t' && nt!='g' && nt!='c') {/*printf ("ERROR not a nucleotide %c\n",c);*/ nt=c;}
	
	return nt;
	}
char genetic_code_hmtDNA (char nt1, char nt2, char nt3)
{
char c;
if ( (nt1!= 'a' && nt1 !='t' && nt1!= 'g' && nt1 !='c') || (nt2!= 'a' && nt2 !='t' && nt2!= 'g' && nt2 !='c') || (nt3!= 'a' && nt3 !='t' && nt3!= 'g' && nt3 !='c') )
		{c='?';
		/*printf("genetic code error ? %c%c%c\n",nt1,nt2,nt3);*/
		}
else if (nt1 == 'a')
		{
		
		if (nt2 == 'c')
		
			{
			c='T';
			}
		else if (nt2 == 'g')
			{
			
			if (nt3 =='t' || nt3 =='c') {c='S';}
			
			else  
				{c='B';
				/*printf ("STOP codon !! \n");*/}
			}
		else if (nt2 == 't')
			{
			if (nt3 =='t' || nt3 =='c') {c='I';}
			else c='M';
			}
		else if (nt2 == 'a')
			{
			if (nt3 =='t' || nt3 =='c') {c='N';}
			else c='K';
			}
		}
		
else if (nt1 == 'c')
	
		{
		
		if (nt2 == 'c')
		
			{
			c='P';
			}
		else if (nt2 == 'g')
			{
			
			c='R';
				
			}
		else if (nt2 == 't')
			{
			c='L';
			}
		else if (nt2 == 'a')
			{
			if (nt3 =='t' || nt3 =='c') {c='H';}
			else c='Q';
			}
		}
else if (nt1 == 'g')
	
		{
		
		if (nt2 == 'c')
		
			{
			c='A';
			}
		else if (nt2 == 'g')
			{
			
			c='G';
				
			}
		else if (nt2 == 't')
			{
			c='V';
			}
		else if (nt2 == 'a')
			{
			if (nt3 =='t' || nt3 =='c') {c='D';}
			else c='E';
			}
		}	
			
else if (nt1 == 't')
	
		{
		
		if (nt2 == 'c')
		
			{
			c='S';
			}
		else if (nt2 == 'g')
			{
			
			if (nt3 =='t' || nt3 =='c') {c='C';}
			else c='W';
				
			}
		else if (nt2 == 't')
			{
			if (nt3 =='t' || nt3 =='c') {c='F';}
			else c='L';
			}
		else if (nt2 == 'a')
			{
			if (nt3 =='t' || nt3 =='c') {c='Y';}
			else { c='B';
					/*printf("STOP codon ! \n");*/}
			}
		}	
		
else {
	printf("%c%c%c \n",nt1,nt2,nt3);
	exit(1);}
return c;		
}	
double ** order (double *tab, int nb) /*prends des donnees appariees et les classes par ordre croissant en gardant leur rang*/
{
int i,j,z,min,max;
double **rangg;
rangg = ftaballoc (nb+1,2);
if (rangg == NULL) {printf("not enough memory to allocate the rank table !!\n"); exit (1);}
rangg[0][0]=tab[0];
rangg[0][1]=0;/*garde la position initiale de la valeur*/
for (i=0;i<(nb-1);i++)
	{
	/*printf("%i %f\n",i,tab[i]);*/
	
	if ( (tab[i+1]>=rangg[i][0])  || (tab[i+1]<=rangg[0][0])   )
		{
		if (tab[i+1]>=rangg[i][0]) 
			{
			rangg[i+1][0]=tab[i+1];
			rangg[i+1][1]=i+1;
			}
		if (tab[i+1]<=rangg[0][0])
			{
			for (z=i; z>=0; z--)
				{
				rangg[z+1][0]=rangg[z][0];
				rangg[z+1][1]=rangg[z][1];
				}
			rangg[0][0]=tab[i+1];
			rangg[0][1]=i+1;
			}	
		
		}
		
	else
		{			
		j=ceil ((i)/2);
		max=i;
		min=0;
		/*printf("%i\n",j);*/
		while (  (tab[i+1] > rangg[j][0])  || (tab[i+1] <= rangg[j-1][0]) )
			{
			
			if (tab[i+1] > rangg[j][0]) { 
										min=j;
										j= ( (ceil (j+(max-j)/2)) < j+1 ) ? j+1 : (ceil (j+(max-j)/2)) ; 
										/*printf("%i\n",j);*/
										}
			
			if (tab[i+1] <= rangg[j-1][0]) {
										max=j-1;
										 j= ( (floor (j-(j-min)/2)) >j-1 ) ? j-1 : (floor (j-(j-min)/2)); 
										 /*printf("%i\n",j);*/
										 }
			
			}
			
		/* rank[j-1][0] < tab[i+1] <= rank[j][0] */
		for (z=i; z>=j; z--)
				{
				rangg[z+1][0]=rangg[z][0];
				rangg[z+1][1]=rangg[z][1];
				}
		rangg[j][0]=tab[i+1];
		rangg[j][1]=i+1;
		
		} 
	
	
	
	}
return rangg;
}
int ** order_int (int *tab, int nb)
{
int i,j,z,min,max;
int **rangg;
rangg = itaballoc (nb,2);
if (rangg == NULL) {printf("\nnot enough memory to allocate the rank table !!\n"); exit (1);}
rangg[0][0]=tab[0];
rangg[0][1]=0;
for (i=0;i<(nb-1);i++)
	{
	/*printf("%i %f\n",i,tab[i]);*/
	
	if ( (tab[i+1]>=rangg[i][0])  || (tab[i+1]<=rangg[0][0])   )
		{
		if (tab[i+1]>=rangg[i][0]) 
			{
			rangg[i+1][0]=tab[i+1];
			rangg[i+1][1]=i+1;
			}
		if (tab[i+1]<=rangg[0][0])
			{
			for (z=i; z>=0; z--)
				{
				rangg[z+1][0]=rangg[z][0];
				rangg[z+1][1]=rangg[z][1];
				}
			rangg[0][0]=tab[i+1];
			rangg[0][1]=i+1;
			}	
		
		}
		
	else
		{
		j=ceil ((i)/2);
		max=i;
		min=0;
		/*printf("%i\n",j);*/
		while (  (tab[i+1] > rangg[j][0])  || (tab[i+1] <= rangg[j-1][0]) )
			{
			
			if (tab[i+1] > rangg[j][0]) { 
										min=j;
										j= ( (ceil (j+(max-j)/2)) < j+1 ) ? j+1 : (ceil (j+(max-j)/2)) ; 
										/*printf("%i\n",j);*/
										}
			
			if (tab[i+1] <= rangg[j-1][0]) {
										max=j-1;
										 j= ( (floor (j-(j-min)/2)) >j-1 ) ? j-1 : (floor (j-(j-min)/2)); 
										 /*printf("%i\n",j);*/
										 }
			
			}
			
		/* rank[j-1][0] < tab[i+1] <= rank[j][0] */
		for (z=i; z>=j; z--)
				{
				rangg[z+1][0]=rangg[z][0];
				rangg[z+1][1]=rangg[z][1];
				}
		rangg[j][0]=tab[i+1];
		rangg[j][1]=i+1;
		
		} 
	
	
	
	}
return rangg;
}
/*returns the ranks of two ordered coupled arrays from order and order_int*/
double ** rank( double **mat_1, int **mat_2, int nb) /* mat1[][0] contains the ordered LD values mat1[][1] the nb of the measure */
												/*mat2[][0] the ordered values of distances mat2[][1] the nb of the measures*/
{
int i,j,sum1,sum2,repet1,repet2;
double **rg;
i=nb;
if ( (rg=ftaballoc(i,2)) == NULL) {printf("\nallocation pb of rank matrice !\n"); exit(1);}
repet1=0;
sum1=0;
repet2=0;
sum2=0;
for (i=0;i<(nb-1);i++)
	{
	
	/*calcul des rangs de mat1*/
	if ( mat_1[i][0]!=mat_1[i+1][0]  &&  mat_1[i][0]<0.99999)
		{
		
		if (repet1 ==0)
			{
			rg[(int) (mat_1[i][1])][0]=i; /*le rang de la ieme valeur est range dans la ieme case*/
			/*printf("rg[%d] %.1f %i\n",(int) (mat1[i][1]),rg[(int) (mat1[i][1])][0],i);*/
			}
		else
			{
			++repet1;
			sum1+=i;
			for (j=i; j>(i-repet1); j--)
				{
				rg[(int) (mat_1[j][1])][0]= (double) (sum1)/ (double) (repet1);
				/*printf("rg[%d] %.1f %i\n",(int) (mat1[j][1]),rg[(int) (mat1[j][1])][0],j);*/
				}
			
			repet1=0;
			sum1=0;
			}
		}
	else 
		{++repet1;
		sum1+=i;
		}
	
	/*calcul des rangs de mat2*/
	if (mat_2[i][0] != mat_2[i+1][0])
		{
		if (repet2 ==0)
			{rg[mat_2[i][1]][1]=i;}
		else 
			{
			++repet2;
			sum2+=i;
			for (j=i; j>(i-repet2); j--)
				{
				rg[mat_2[j][1]][1]=(double) (sum2)/(double) (repet2);
				}
			
			repet2=0;
			sum2=0;
			}
		}
	else 
		{++repet2;
		sum2+=i;
		}
	}
/*for the last line*/
	
if (repet1 ==0)
			{rg[(int) (mat_1[nb-1][1])][0]=nb-1;}
else 
			{
			++repet1;
			sum1+=nb-1;
			for (j=nb-1; j>(nb-1-repet1); j--)
				{
				rg[(int) (mat_1[j][1])][0]=(double) (sum1)/(double) (repet1);
				}
			
			repet1=0;
			sum1=0;
			}
			
if (repet2 ==0)
			{rg[(int) (mat_2[nb-1][1])][1]=nb-1;}
else 
			{
			++repet2;
			sum2=nb-1;
			for (j=nb-1; j> (nb-1-repet2); j--)
				{
				rg[(int) (mat_2[j][1])][1]=(double) (sum2)/(double) (repet2);
				}
			
			repet2=0;
			sum2=0;
			}
return rg;
}
double spearman (double ** rang, int nb)
{
int i;
double rho;
long double Srisi,Sri2,Ssi2;
Srisi=0;
Sri2=0;
Ssi2=0;
/*printf("%f %f %f\n", Srisi,Sri2,Ssi2);*/
for (i=0;i<nb;i++)
	{
	Srisi=Srisi+(rang[i][0]+1)*(rang[i][1]+1)-(nb+1)*(nb+1)/4;
	Sri2=Sri2+(rang[i][0]+1)*(rang[i][0]+1)-(nb+1)*(nb+1)/4;
	Ssi2=Ssi2+(rang[i][1]+1)*(rang[i][1]+1)-(nb+1)*(nb+1)/4;
	}
rho=Srisi/(sqrt(Sri2*Ssi2));
	
return rho;	
}
double Pearson (double *LD, int *dist, int nb)
{
int i;
double rho_p;
double Sxiyi,Sxi2,Syi2,x_mean,y_mean;
Sxiyi=0;
Sxi2=0;
Syi2=0;
x_mean=y_mean=0;
/*printf("%f %f %f\n", Srisi,Sri2,Ssi2);*/
for (i=0;i<nb;i++)
	{
	x_mean+=LD[i];
	y_mean+=dist[i];
	}
x_mean=x_mean/nb;
y_mean=y_mean/nb;
for (i=0;i<nb;i++)
	{
	Sxiyi=Sxiyi+(LD[i]-x_mean)*(dist[i]-y_mean);
	Sxi2=Sxi2+(LD[i]-x_mean)*(LD[i]-x_mean);
	Syi2=Syi2+(dist[i]-y_mean)*(dist[i]-y_mean);
	}
rho_p=Sxiyi/(sqrt(Sxi2)*sqrt(Syi2));
	
return rho_p;	
}
/*
double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
*/

double permute_Pearson (double *LD, int *pos, int permut, double rho, int nb)
{
int i,j,**p, **newpos,couple, **newdist, *newdist2,k,nbcase,rag,compteur,max,min;
double gg, P_value,P_value2, *P, **tmp, *R,zz;
//long idum;

P_value=1.0;

couple=nb*(nb-1)/2;
p=itaballoc(nb,2);
	if(p==NULL) {printf("\nnot enough memory to permute P\n");exit(1);}
newpos=itaballoc(nb,2);
	if(newpos==NULL) {printf("\nnot enough memory to permute P\n");exit(1);}
	if( (newdist=itaballoc(couple+1,2)) ==NULL) {printf("\nnot enough memory to permute P\n");exit(1);}
newdist2=calloc(couple+1,sizeof(int));
	if(newdist2==NULL) {printf("\nnot enough memory to permute P\n");exit(1);}
R=calloc(permut+1, sizeof(double));
	if(R==NULL) {printf("\nnot enough memory to permute P\n");exit(1);}
tmp=ftaballoc(permut,2);
	if(tmp==NULL) {printf("\nnot enough memory to permute P\n");exit(1);}

for (i=0;i<permut;i++)
	{
	nbcase=nb;
	for (j=0;j<nb;j++)
		{
		p[j][0]=pos[j];
		p[j][1]=j; /*garde en memoire la position initiale*/
		}
	
	j=0;
	while (nbcase!=0)
		{
		//gg=ran1(&idum);
		//rag=gg*nbcase;
		rag=Rand(nbcase);
		newpos[j][0]=p[rag][0]; /*les nlles positions sont choisies a partir du tableau p*/

		newpos[j][1]=p[rag][1]; /*memoire de la position initiale, pour retrouver le ld associe !
		printf("%d\t %d\n",newpos[j][0],newpos[j][1]);*/
		
		for (k=rag;k<(nbcase-1);k++)
			{
			p[k][0]=p[k+1][0];
			p[k][1]=p[k+1][1];
			}
		nbcase--;
		j++;
		}
		
	/*calcul des nvlles distances */
	
	compteur=0;
	for(j=0;j<nb-1;j++)
		{
		for(k=j+1;k<nb;k++)
			{
			newdist[compteur][0]=( abs(newpos[j][0]-newpos[k][0]) < (16569-abs(newpos[k][0]-newpos[j][0])) ) ? abs(newpos[k][0]-newpos[j][0]) : (16569-abs(newpos[k][0]-newpos[j][0]));
			newdist2[compteur]=newdist[compteur][0];
			
			
			compteur++;
			}
		}
	/*ordonner le nouveau vecteur distance*/
	
	R[i]=Pearson(LD,newdist2,couple);
	/*printf("R[%i] is %f \n",i,R[i]);*/
	}
tmp=order(R, permut);
zz=0;


for(i=0;i<permut;i++)
	{
	if (rho<tmp[i][0])
		{
		
			//printf("i=%d rho=%f R[i]=%f\n",i,rho,R[i]);
			P_value=(float)(i)/ (float) (permut);
			//printf("\n*P_value %f\n",P_value);
			i=permut;
			
		}
	

	}

return P_value;
}
double permute_spearman2 (double ** LDordonne, int *pos, int permut, double rho, int nb) /*permute les positions ordonne et calcule Prob(rho<rhoobs)*/
{
int i,j,nbcase=nb, rag,k, compteur, min, max,ff;
int **p, **newpos,couple, **newdist, *newdist2, Zpos, **mat22,zz;
double gg, P_value, P_value2,*R, **new_rank, **tmp;
//long idum;
couple=nb*(nb-1)/2;
p=itaballoc(nb,2);
	while (p==NULL) {printf("\nout of memory permute spearman\n"); exit(1);}
newpos=itaballoc(nb,2);
	while (newpos==NULL) {printf("\nout of memory permute spearman\n"); exit(1);}
newdist=itaballoc(couple,2);
	while (newpos==NULL) {printf("\nout of memory permute spearman\n"); exit(1);}
	
newdist2=calloc(couple,sizeof(int));
	while (newpos==NULL) {printf("\nout of memory permute spearman\n"); exit(1);}
	
mat22=itaballoc(couple,2);
	while (newpos==NULL) {printf("\nout of memory permute spearman\n"); exit(1);}
if (p==NULL || newpos==NULL || newdist==NULL || newdist2==NULL || mat22==NULL)
		{
		printf("\nnot enough memory to allocate integer tables in permute_spearman2\n");
		exit(1);
		}
		
R=calloc(permut, sizeof(double));
new_rank=ftaballoc(couple,2);
tmp=ftaballoc(permut,2);
if (R==NULL || new_rank==NULL || tmp==NULL )
	{
		printf("\nnot enough memory to allocate float tables in permute_spearman2\n");
		exit(1);
		}
for (i=0;i<permut;i++)
	{
	
	/*if (i%(permut/5)==0) printf("\t%i",i);*/
	
	nbcase=nb;
	
	for (j=0;j<nb;j++)
	{
	p[j][0]=pos[j];
	p[j][1]=j; /*garde en memoire la position initiale*/
	}
	
	j=0;
	while (nbcase!=0)
		{
		//gg=ran1(&idum);
		//rag=gg*nbcase;
		rag=Rand(nbcase);
		newpos[j][0]=p[rag][0]; /*les nlles positions sont choisies a partir du tableau p*/
		
		newpos[j][1]=p[rag][1]; /*memoire de la position initiale, pour retrouver le ld associe !
		printf("%d\t %d\n",newpos[j][0],newpos[j][1]);*/
		
		for (k=rag;k<(nbcase-1);k++)
			{
			p[k][0]=p[k+1][0];
			p[k][1]=p[k+1][1];
			}
		nbcase--;
		j++;
		}
		
	/*calcul des nvlles distances */
	
	compteur=0;
	for(j=0;j<nb-1;j++)
		{
		for(k=j+1;k<nb;k++)
			{
			newdist[compteur][0]=( abs(newpos[j][0]-newpos[k][0]) < (16569-abs(newpos[k][0]-newpos[j][0])) ) ? abs(newpos[k][0]-newpos[j][0]) : (16569-abs(newpos[k][0]-newpos[j][0]));
			newdist2[compteur]=newdist[compteur][0];
			/*printf("j %i k %i %i %i %i\n",j,k,newdist[compteur][0],newpos[j][0],newpos[k][0]);*/
			
			min= (newpos[j][1] < newpos[k][1]) ? newpos[j][1] : newpos[k][1];
			max= (newpos[j][1] > newpos[k][1]) ? newpos[j][1] : newpos[k][1];
			
			Zpos=0;

			Zpos=Zpos+max-min-1;

			newdist[compteur][1]=Zpos;
			
			compteur++;
			}
		}
	/*ordonner le nouveau vecteur distance*/
	
	mat22=order_int(newdist2,couple);
	
	new_rank=rank( LDordonne, mat22, couple);
	
	R[i]=spearman(new_rank,couple);

	}
tmp=order(R, permut);

P_value=0;
zz=0;

for(i=0;i<permut;i++)
	{
	 if (rho<tmp[i][0])
		{
		
			/*printf("%d %f %f\n",i,rho,R[i]);*/
		P_value=(float)(i)/ (float) (permut);
		i=permut;
			
		}
		
	if (zz==0 && rho>0 && rho>tmp[permut-i-1][0])
		{
		
			/*printf("%d %f %f\n",i,rho,R[i]);i=permut;*/
			P_value2=(float)(i+1)/ (float) (permut);
			printf(" P2 %f\t",P_value2);
			zz=1;
		}
	}
		
printf("P_value   =%f\t",P_value);
free(R);
for (i=0;i<nb;i++)
	{
	free(p[i]);
	free(newpos[i]);
	}
	
for (i=0;i<couple;i++)
	{
	free(newdist[i]);
	free(mat22);
	free(new_rank);
	free(tmp);
	}
return P_value;
}
/* permutes permut times the rank of the distances and calculates rhos, orders rho and gives the percent of rhos smaller than the observed rho*/
int is_nt (char nt)
{
int c;
c=0;
if (tolower(nt)=='a' || tolower(nt)=='t' || tolower(nt)=='g' || tolower(nt)=='c')
	{c=1;}
else c=0;
return c;
}

void SetSeed (int seed)
{
   z_rndu = 170*(seed%178) + 137;
   w_rndu=seed;
}
#ifdef FAST_RANDOM_NUMBER
double rndu (void)
{
   w_rndu *= 127773;
   return ldexp((double)w_rndu, -32);
}
#else
double rndu (void)
{
   static int x_rndu=11, y_rndu=23;
   double r;
   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
   return (r-(int)r);
}
#endif
/********************** Rand ****************************/
/* it returns number in range 0...max-1                   */
static int Rand(int max)
{
	double	rd;
/*	rd = rand()/32767.0; /* returns a number in range 0.0-1.0 */
	rd = rndu();
	return (floor(rd*max));
}
/********************* PrintDate ***********************/
static void PrintDate ()
{
	time_t now;
	char *date;

	#ifdef macintosh
		printf("Macintosh OS\n");
	#endif

	now=time(NULL);
	date= ctime(&now);
	printf("%s",date);
}

