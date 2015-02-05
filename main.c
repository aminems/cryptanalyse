/*
Just a base to do some experiments 
Mohamed Amine Bergach

inspired from : Communication Theory of Secrecy Systems
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LG_MESS_MAX 300

//declaration des variables globales

int nbrAlphabe=0;
int lenghtOfWord=0;
char* alphabet=NULL;
double* probAlpha=NULL;
char **tabClefs;
char* combi;
int indiceTabClefs=0;
int indiceTabCombi=0;
double **pYsK;
double **pKsY;
double **pYX;
double *pC;
double *pYeC;
double info=0;
double nLatta;
char mess[LG_MESS_MAX];
char messCrypt[LG_MESS_MAX];
char messDecrypt[LG_MESS_MAX];
double **tabAtta;
double *produitProb;
int *indice;
//factoriel
double fact(int n)
{
	double f = 1;
	while (n > 1) f *= n--;
	return f;
}

void probaYsK()
{
	int t,g,lk,i,l,j;
	int lenY=(int)pow(nbrAlphabe, lenghtOfWord)*lenghtOfWord;
	int N=(int)fact(nbrAlphabe);
	int Nl=(int)pow(nbrAlphabe, lenghtOfWord);
	pYsK=(double **) malloc(N * sizeof(double*));
	for( t = 0; t < N; t++)
		pYsK[t] =(double *) malloc(Nl * sizeof(double));
	for ( g=0;g<Nl;g++) 
	{
		for ( lk=0; lk<fact(nbrAlphabe); lk++)
		{
			for ( i=0;i<lenY;i=i+lenghtOfWord)
			{
				for ( l=i; l<lenghtOfWord; l++) 
				{
					for ( j=0;j<nbrAlphabe;j++)
					{
						if(tabClefs[lk][j]==combi[i])
						{
							if(pYsK[lk][g]==0)
								pYsK[lk][g]=probAlpha[j];
							else
							pYsK[lk][g]*=probAlpha[j];
						}
					}
				}
			}
		}
	}
}


void probaKsY(int nbC, int nbL)
{
	int t,i,j;
	int N=(int)fact(nbrAlphabe);
	int Nl=(int)pow(nbrAlphabe, lenghtOfWord);
	pKsY=(double **) malloc(N * sizeof(double*));
	for( t = 0; t < N; t++)
		pKsY[t] =(double *) malloc(Nl * sizeof(double));
	for ( i=0;i<=(fact(nbL)-1);i++)
	{
		for ( j=0;j<=(nbC-1);j++)
		{
			pKsY[i][j]=pYsK[i][j]/fact(nbL-1);
		}
	}
}

void probaYX(int LW, int nC, int nL)
{
	int a, b, c, d, e, i, j, k, t;
	int comb3LettersEqual[nC];
	int comb2FirstLettersEqual[nC];
	int combFirstAndLastLettersEqual[nC];
	int comb2LastLettersEqual[nC];
	int combAllDifferentLetters[nC]; 
	int combLettersEqual[nC];
	int combLettersDifferent[nC];
	int Indice;
	for (Indice = 0; Indice < nC; Indice++) {
		comb3LettersEqual[Indice] = 0;
		comb2FirstLettersEqual[Indice] = 0;
		combFirstAndLastLettersEqual[Indice] = 0;
		comb2LastLettersEqual[Indice] = 0;
		combAllDifferentLetters[Indice] = 0;
		combLettersEqual[Indice] = 0;
		combLettersDifferent[Indice] = 0;
	} 
	
	pYX=(double **) malloc(nC * sizeof(double*));
	for( t = 0; t < nC; t++)
		pYX[t] =(double *) malloc(nC * sizeof(double));
	
	if (LW==3)
	{ 
		a=b=c=d=e=0;
		for (k=0;k<=nC-1;k++)
		{
			if ((combi[k]==combi[k+1]) &&  (combi[k]==combi[k+2]))
			{
				a+=1;
				comb3LettersEqual[a-1]=k;
			}
			else if ((combi[k]==combi[k+1]) &&  (combi[k]!=combi[k+2]))
			{
				b+=1;
				comb2FirstLettersEqual[b-1]=k;
			}
			else if ((combi[k]==combi[k+2]) &&  (combi[k]!=combi[k+1]))
			{
				c+=1;
				combFirstAndLastLettersEqual[c-1]=k;
			}
			else if ((combi[k+1]==combi[k+2]) &&  (combi[k+2]!=combi[k+1]))
			{
				d+=1;
				comb2LastLettersEqual[d-1]=k;
			}
			else if (((combi[k+1]!=combi[k+2]) &&  (combi[k+2]!=combi[k])) && (combi[k]!=combi[k+1]))
			{
				e+=1;
				combAllDifferentLetters[e-1]=k;
			}
		}
		for (i=0;i<=a-1;i++)
		{
			for (j=0;j<=a-1;j++)
			{
				pYX[comb3LettersEqual[i]][comb3LettersEqual[j]]=1/fact(nL);

			}
		}
		for (i=0;i<=b-1;i++)
		{
			for (j=0;j<=b-1;j++)
			{
				pYX[comb2FirstLettersEqual[i]][comb2FirstLettersEqual[j]]=1/fact(nL);
			}
		}
		for (i=0;i<=c-1;i++)
		{
			for (j=0;j<=c-1;j++)
			{
				pYX[combFirstAndLastLettersEqual[i]][combFirstAndLastLettersEqual[j]]=1/fact(nL);
			}
		}
		for (i=0;i<=d-1;i++)
		{
			for (j=0;j<=d-1;j++)
			{
				pYX[comb2LastLettersEqual[i]][comb2LastLettersEqual[j]]=1/fact(nL);
			}
		}
		for (i=0;i<=e-1;i++)
		{
			for (j=0;j<=e-1;j++)
			{
				pYX[combAllDifferentLetters[i]][combAllDifferentLetters[j]]=1/fact(nL);
			}
		}
	}
	else if (LW==2)
	{ 
		a=b=0;
		for (k=0;k<=nC-1;k++)
		{
			if (combi[k]==combi[k+1])
			{
				a+=1;
				combLettersEqual[a-1]=k;
			}
			else if (combi[k]!=combi[k+1])
			{
				b+=1;
				combLettersDifferent[b-1]=k;
			}
		}
		for (i=0;i<=a-1;i++)
		{
			for (j=0;j<=a-1;j++)
			{
				pYX[combLettersEqual[i]][combLettersEqual[j]]=1/fact(nL);
			}
		}
		for (i=0;i<=b-1;i++)
		{
			for (j=0;j<=b-1;j++)
			{
				pYX[combLettersDifferent[i]][combLettersDifferent[j]]=1/fact(nL);
			}
		}
	}
	else if (LW==1)
	{
		for (i=0;i<=nL-1;i++)
		{
			for (j=0;j<=nL-1;j++)
			{
				pYX[combLettersEqual[i]][combLettersEqual[j]]=1/fact(nL);
			}
		}
	}
}



void probComb(int lW, int nL,int nC)
{ 
	int k,i,j;
	pC=(double *)malloc(((int)pow(nbrAlphabe, lenghtOfWord)) * sizeof(double));
	if(lW!=1)
	{
		for (k=0;k<=(nC-1);k++)
		{ 
			pC[k]=1;
			for (i=lW*k;i<=(((k+1)*lW)-1);i++)
			{
				for (j=0;j<=nL;j++)
				{
					if (combi[i]==alphabet[j])
					{
						pC[k]=pC[k]*probAlpha[j];
					}
				}
			}
		}
	}
	else {
		for (i=0;i<=nL-1;i++)
		{
			pC[i]=probAlpha[i];
		}
	}
}




void probYegalComb(int nC)
{ 
	int i,j;
	pYeC=(double *)malloc(((int)pow(nbrAlphabe, lenghtOfWord)) * sizeof(double));
	pYeC[0]=0;
	for(i=0;i<=nC-1;i++)
	{
		for (j=0;j<=nC-1;j++)
		{
			pYeC[i]=pYeC[i]+pC[j]*pYX[i][j];
		}
	}
}



void infoYK(int nC)
{
	double hY=0;
	double hYsK=0;
	int i=0;
	int j=0;
	for( i=0;i<nC;i++)
	{
		hY+=pYeC[i]*log2(pYeC[i]);
	}
	int N=(int)fact(nbrAlphabe);
	int NN=(int)pow(nbrAlphabe, lenghtOfWord);
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < NN; j++)
		{
			hYsK+=pYsK[i][j]*log2(pYsK[i][j]);
		}
	}
	info=hY-hYsK;
}

void nbLettersForDecipher()
{ 
	nLatta=log2(fact(nbrAlphabe))/info;
	/*nLatta=round(nLatta);
	double rest=(int)nLatta%lenghtOfWord;
	nLatta+=lenghtOfWord-rest;*/
}





void generateMessage()
{
	int k=0;
	srand(nbrAlphabe);
	while (k<LG_MESS_MAX)
	{
		mess[k]=alphabet[(int)(rand()%nbrAlphabe)];
		k+=1;
	}
}

void crypt(int keyNb) 
{
	int i;
	for (i=0;i<=LG_MESS_MAX;i++)
	{ 
		if (mess[i]==alphabet[0])
		{ 
			messCrypt[i]=tabClefs[keyNb][0];}
		else if (mess[i]==alphabet[1])
		{ 
			messCrypt[i]=tabClefs[keyNb][1];}
		else if (mess[i]==alphabet[2])
		{ 
			messCrypt[i]=tabClefs[keyNb][2];}
		else if (mess[i]==alphabet[3])
		{ 
			messCrypt[i]=tabClefs[keyNb][3];}
	}
}

void decrypt(int keyNb)
{
	int i,j;
	for (i=0;i<=LG_MESS_MAX;i++)
	{ 
		for (j=0;j<=(nbrAlphabe-1);j++)
		{
			if (tabClefs[keyNb][j]==messCrypt[i])
				messDecrypt[i]=alphabet[j];
		}
	}
}

int chooseRandomlyAKey()
{
	int keynb;  
	srand(nbrAlphabe);
	keynb=(int)rand()%(int)fact(nbrAlphabe);
	return keynb;
}

void cryptingANDdecrypting(int choice, int keyNumber)
{
	int userkey, userchoice, counter=0;
	
	if (choice==0)
	{ 
		while (counter==0)
		{
			printf("The message given will be encrypted.\n If you want to choose the key for the encryption enter 0, if not enter 1 and a key will be pseudo-randomly chosen.\n Choice :\n"); 
			scanf("%d", &userchoice);
			if (userchoice==0)
			{ 
				printf("Please enter the number of the key you want to use for the encryption: \n");
				scanf("%d", &userkey);
				crypt(userkey);
				counter=1;
			}
			else if (userchoice==1)
			{
				userkey=chooseRandomlyAKey();
				crypt(userkey);
				counter=1;
			}
			else printf("Error! \n The number given is neither 0 nor 1. Please enter the correct value. \n");
		}
	}
	else if (choice==1)
	{
		printf("The message given will be decrypted with the key %d . \n", keyNumber);
		decrypt(keyNumber);
	}
}

void tabAttack()
{
	int t, k, j;
	int i=0;
	int g=0;
	int N=(int)fact(nbrAlphabe);
	int Nl=(int)nLatta;
	tabAtta=(double **) malloc(N * sizeof(double*));
	for( t=0; t<N; t++)
		tabAtta[t] =(double *) malloc(Nl * sizeof(double));
	while(i<nLatta)
	{
		for ( k=0; k<=(pow(nbrAlphabe, lenghtOfWord)-1);k++)
		{
			for (i=lenghtOfWord*k; i<(((k+1)*lenghtOfWord)); i++)
  		 	{
				if (combi[i]==messCrypt[i])
 			 	{
					g+=1;
				}
			}   
			if (g==lenghtOfWord)
		 	{
				for(  j=0; j<=fact(nbrAlphabe)-1; j++)
				{
					tabAtta[j][k]=pKsY[j][k];
				}
			}
			else g=0;
		}
	}
}

void prodProb()
{
	int i,j;
	produitProb=(double *)malloc(((int)fact(nbrAlphabe) * sizeof(double)));
	produitProb[0]=1;
	for ( i=0; i<=fact(nbrAlphabe)-1; i++)
	{
		for ( j=0; j<=nLatta-1 ; j++)
		{
			produitProb[i]=produitProb[i]*tabAtta[i][j];
		}
	}
} 

void getTheGoodKey()
{
	int i;
	for ( i=0;i<=(fact(nbrAlphabe)-1);i++)
	{
		if (produitProb[indice[0]]<produitProb[i])
		{
			indice[0]=i;
		}
	}
	for (i=0;i<=(fact(nbrAlphabe)-1);i++)
	{
		if ((produitProb[indice[0]]==produitProb[i]) && (indice[0]!=i))
		{
			indice[1]=i;
			printf("Existence of two keys which were possibly used to crypt:\n key k%d and key k%d .\n", indice[0], indice[1]);
		}
	}
}



//génére les combinaison
void combinaison(char sSeq[], int len, char sOut[], int level, int depth)
{
	int i,j;
	for ( i=0; i<len; i++)
	{
		sOut[level-1]= sSeq[i];
		if ( level == depth)
		{
			for ( j=0; j<depth; j++) {
				combi[indiceTabCombi]=sOut[j];
				indiceTabCombi++;
			}
			//printf("%s\n", sOut);
		}
		else
			combinaison(sSeq, len, sOut, level+1, depth);
	}
}

//génére les clefs
void permute(char *string_start, char *p) {
	
	if (*(p+1) == 0) {
		strcpy(tabClefs[indiceTabClefs],string_start);
		//printf("%d %s\n",indiceTabClefs, tabClefs[indiceTabClefs]);
		indiceTabClefs++;
		
	}
	else {
		char *swap;
		for(swap = p; *swap; ++swap) {
			char *same;
			for(same = p; *same != *swap; ++same) {
			}
			if (same == swap) {
				char tmp = *swap;
				*swap = *p;
				*p = tmp;
				permute(string_start, p+1);
				*p = *swap;
				*swap = tmp;
			}
		}
	}
}






int main (int argc, const char * argv[])
{    
	int i=0;
	int N=0;
	double totalProb=0;
	printf("Please enter the number of letters of the alphabet:\n");
	scanf("%d",&nbrAlphabe);
	
	alphabet=(char *)malloc((nbrAlphabe+1) * sizeof(char));
	probAlpha=(double *)malloc((nbrAlphabe+1) * sizeof(double));
	N=(int)fact(nbrAlphabe);
	tabClefs=(char **) malloc(N * sizeof(char*));
	if(tabClefs == NULL){
		free(tabClefs);
		printf("Memory allocation failed while allocating for dim[].\n");
		exit(-1);
	}
	for(i = 0; i < N; i++)
	{
		tabClefs[i] =(char *) malloc((nbrAlphabe+1) * sizeof(char));
		if(NULL == tabClefs[i])
		{
			free(tabClefs[i]);
			printf("Memory allocation failed while allocating for dim[x][].\n");
			exit(-1);
		}
	}
	while (totalProb!=1){
		for(i=0;i<nbrAlphabe;i++)
		{
			printf("\n********Please enter a letter:\n");
			scanf("\n%c",&alphabet[i]);
			printf("\n********Please enter the probability of  %c :\n",alphabet[i]);
			scanf("\n%lf",&probAlpha[i]);
			while(probAlpha[i]>=1 || probAlpha[i]<=0) {
				printf("\nERROR PROBA Please REenter the probability of  %c :\n",alphabet[i]);
				scanf("\n%lf",&probAlpha[i]);
			}
			totalProb+=probAlpha[i];
		}
		if (totalProb!=1) {
			printf("\nERROR PROBA Please REenter ALL the probabilities:\n");
			totalProb=0;
		}
	}
	
	for(i=0;i<nbrAlphabe;i++)
	{
		printf("probability of: %c %lf\n",alphabet[i],probAlpha[i]);
	}
	printf("Please enter the length of words:\n");
	scanf("%d",&lenghtOfWord);     
	printf("\n");
	
	int NombreDeCombinaisons=(int)pow(nbrAlphabe, lenghtOfWord);
	combi=(char *)malloc(((int)pow(nbrAlphabe, lenghtOfWord)*lenghtOfWord) * sizeof(char));
	permute(alphabet, alphabet);
	printf("ALL possible KEYS : \n");
	for(int f = 0; f < N; f++)
	{
		for(int ff = 0; ff < nbrAlphabe; ff++)
		{
			printf("%c",tabClefs[f][ff]);
		}
		printf("  (K%d)\n",f);
	}
	char sOut[lenghtOfWord];
	combinaison(alphabet, nbrAlphabe, sOut, 1, lenghtOfWord);
	printf("ALL possible Combinations : \n");
	for(int f = 0; f <NombreDeCombinaisons*lenghtOfWord; f=f+lenghtOfWord)
	{
		for(int ff = f; ff < lenghtOfWord+f; ff++)
		{
			printf("%c",combi[ff]);
		}
		printf("\n");
	}
	probaYsK();
	printf("probability of Y/K : \n");
	int NN=(int)pow(nbrAlphabe, lenghtOfWord);
	for(int f = 0; f <N; f++)
	{
		for(int ff = 0; ff < NN; ff++)
		{
			printf(" %lf ",pYsK[f][ff]);
		}
		printf("\n");
	}
	probaKsY(NombreDeCombinaisons, nbrAlphabe);
	printf("probability of K/Y : \n");
	for(int f = 0; f <N; f++)
	{
		for(int ff = 0; ff < NN; ff++)
		{
			printf(" %lf ",pKsY[f][ff]);
		}
		printf("\n");
	}
	probaYX(lenghtOfWord,NombreDeCombinaisons,nbrAlphabe);
	printf("probability of YX : \n");
	for(int f = 0; f <NombreDeCombinaisons; f++)
	{
		for(int ff = 0; ff < NombreDeCombinaisons; ff++)
		{
			printf(" %lf ",pYX[f][ff]);
		}
		printf("\n");
	}
	probComb(lenghtOfWord,nbrAlphabe,NombreDeCombinaisons);
    printf("probability of Combinations : \n");
	for(int f = 0; f <NombreDeCombinaisons; f++)
	{
			printf(" %lf ",pC[f]);
	}
	printf("\n");
	probYegalComb(NombreDeCombinaisons);
	printf("probability of Y=Comb : \n");
	for(int f = 0; f <NombreDeCombinaisons; f++)
	{
		printf(" %lf ",pYeC[f]);
	}
	printf("\n");
	infoYK(NombreDeCombinaisons);
	printf("Mutual information amount:  %lf \n",info);
	nbLettersForDecipher();
	printf("Number of letter for Decipher:  %lf \n",nLatta);
	generateMessage();
	printf("Original Message: ");
	for(int f = 0; f <LG_MESS_MAX; f++)
	{
		printf("%c",mess[f]);
	}
	printf("\n");
	cryptingANDdecrypting(0,0);
	printf("Message to Decipher: ");
	for(int f = 0; f <LG_MESS_MAX; f++)
	{
		printf("%c",messCrypt[f]);
	}
	printf("\n");
	tabAttack();
	printf("Attack Table: \n");
	int Nla=(int)nLatta;
	for(int f = 0; f <N; f++)
	{
		for(int ff = 0; ff < Nla; ff++)
		{
			printf(" %lf ",tabAtta[f][ff]);
		}
		printf("\n");
	}
	prodProb();
	getTheGoodKey();
	cryptingANDdecrypting(1,indice[0]);
	if (indice[1]!=NULL)
	{
		cryptingANDdecrypting(1,indice[1]);
	}
	
	return 0; 
}