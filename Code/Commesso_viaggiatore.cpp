#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdio.h>

#define IA (16807)     
#define IM (2147483647)   
#define AM (1.0L/IM)

using namespace std;

double ran0(long *idum); // Generatore di numeri casuali di tipo LCG con distribuzine uniforme in [0,1]
long idum;                

int main()
{
	int N = 0; // Numero totale delle città intero e quadrato perfetto
	int L = 0; 
	int Nk = 0; // Numero di Sweeps di Montecarlo
	int cont = 0; // Conta le energie accettate 
	int passo = 0; // Conta il numero di volte con cui viene ripetuto l'algoritmo
	int num = 0; // Numero casuale generato per la configurazione iniziale
	int a = 0; // Numero casuale in [0,N-1] generato per lo switch delle posizioni
	int b = 0; // Numero casuale in [0,N-1] generato per lo switch delle posizioni

	double pacc = 0.0; // Probabilità di accettazione
	double amin = 0.0; // Tasso di accettazione sotto il quale l'agoritmo viene interrotto
	double ak = 1.0; // Tasso di accettazione calcolato
	double delta = 0.0; // Parametro per l'incremento della Beta
	double Bl= 0.0;	// Beta della simulazione precedente
	double Bl1= 0.0; // Beta della simulazione corrente
	double Ek1 = 0.0; // Energia della simulazione corrente 
	double Ej = 0.0; // Somma delle energie accettate
	double Ej2 = 0.0; // Somma del quadrato delle energie accettate
	double medEk = 0.0; // Media delle energie accettate
	double medEk2 = 0.0; // Media del quadrato delle energie accettate
	double varEk = 0.0; // Varianza delle energie
	double c = 0; // Numero casuale in [0,1] generato per il confronto con la probabilità di accettazione
	double d = 0; // Numero casuale in [0,1] generato per il confronto con la probabilità di accettazione
	double Eold = 0; // Energia della configurazione precedente
	double Etrial = 0; // Energia della configurazione di prova 
	double deltaE = 0; // Differenza di energia tra le configurazioni

	// Variabili per il calcolo di deltaE
	double am1 = 0; 
	double ap1 = 0;
	double bm1 = 0;
	double bp1 = 0;
	double xa = 0;
	double xap1 = 0;
	double xam1 = 0;
	double ya = 0;
	double yap1 = 0;
	double yam1 = 0;
	double xb = 0; 
	double xbp1 = 0;
	double xbm1 = 0;
	double yb = 0;
	double ybp1 = 0;
	double ybm1 = 0;
	double xxa = 0;
	double xxap1 = 0;
	double xxam1 = 0;
	double yya = 0;
	double yyap1 = 0;
	double yyam1 = 0;
	double xxb = 0; 
	double xxbp1 = 0;
	double xxbm1 = 0;
	double yyb = 0;
	double yybp1 = 0;
	double yybm1 = 0;

	cout << "PROGRAMMA PER LA RISOLUZIONE DEL PROBLEMA DEL COMMESSO VIAGGIATORE" << endl;
	cout << "---------------------------------------------------------------------" << endl;	
	cout << "Inserire il lato del quadrato della griglia: ";
	cin >> L;
	cout << "Inserire il numero di Sweeps di Montecarlo: ";
	cin >> Nk;
	cout << "Inserire il valore del parametro Delta per per l'incremento della Beta: ";
	cin >> delta;
	cout << "Inserire il valore del tasso di accettazione minimo: ";
	cin >> amin;
	cout << "Inserire il valore del seme (intero) per la generazione dei numeri casuali: ";
	cin >> idum;
	cout << "---------------------------------------------------------------------" << endl;

	ofstream grafen;
	grafen.open("Passi_Energia.txt");
	ofstream periniz;
	periniz.open("Percorso_iniziale.txt");
	ofstream perfin;
	perfin.open("Percorso_finale.txt");

	N = L*L;

	vector <int> olddisp (N,0); // Configurazione delle posizioni delle città
	vector <int> newdisp (N,0); // Configurazione delle posizioni delle città
	vector <int> confiniziale (N,0); // Configurazione iniziale delle posizioni delle città 
	vector <int> x (N); // Ascissa della posizione della città nella griglia
	vector <int> y (N); // Ordinata della posizione della città nella griglia
	vector <int> xstart (N); // Ascissa della posizione della città nella griglia
	vector <int> ystart (N); // Ordinata della posizione della città nella griglia

// Genero la configurazione iniziale in modo casuale	
	int i = 0;
	while (i < N)
		{		
 		num = floor(N*(ran0(&idum))); // Genero un numero casuale 
		
		if (olddisp.at(num)==0) // Se la posizione del vettore relativa a num è vuota inserisco il valore i+1
			{
			olddisp.at(num)=i+1;
			i++;
			}
		}
	
	confiniziale = olddisp; // Salvo la configurazione iniziale

	for (int i = 0; i < N ; i++) // Calcolo le posizioni delle città nella griglia
		{
		if ((confiniziale.at(i))%L == 0)
			{
			x.at(i) = L ;
			y.at(i) = floor(confiniziale.at(i)/L);
			}
		else
			{
			x.at(i) = confiniziale.at(i)%L ;
			y.at(i) = floor(confiniziale.at(i)/L) + 1;
			}
		}

	xstart = x;
	ystart = y;

// Ogni ciclo abbasso la temperatura fino a che il tasso di accettazione risulta minore del parametro di riferimento inserito	
	while (ak>amin)
		{
		if (passo == 0) // Per il primo set di simulazioni la temperatura è infinita 
			{Bl=0;}
		else // Per i passi successivi calcolo la varianza delle energie e vario Beta in funzione di questa
			{					
			medEk = Ej/cont;		
			medEk2 = Ej2/cont;
			varEk = medEk2 - medEk*medEk; // Calcolo la varianza
			Ej = 0;
			Ej2 = 0;
			
			if (passo == 1) // Per il secondo passo tramite questa formula
				{Bl = 1/sqrt(varEk);}
			else
				{
				Bl1 = Bl + Bl*(log(1+delta))/(3*Bl*(sqrt(varEk))); // Per i successivi tramite questa formula
				Bl = Bl1;
				}
			}
		cont = 0; // Azzero il contatore delle energie accettate prima di riiziale il giro

	//Simulazione di Monte Carlo per N*Nk volte alla stessa Beta=Bl 			
		for (int k = 0; k < Nk*N; k++)
			{
			newdisp = olddisp; //Copio la vecchia configurazione nella nuova		

 			a = floor(N*(ran0(&idum))); // Genero due numeri casuali interi in [0,N-1]				
			d = ran0(&idum);

			if (d < (1-1/(exp(deltaE*Bl))))
				{b = floor(N*(ran0(&idum)));}
			else
				{
				if (a==N-1) {b=0;}
				else {b=a+1;}
				}

			newdisp.at(a) = olddisp.at(b); // Scambio le posizioni dell'interno del vettore
			newdisp.at(b) = olddisp.at(a);

			// Vecchia Confgurazione
			//______________________________________________________________
			// Posizione a
			if ((olddisp.at(a))%L == 0)
				{
				xxa = L ;
				yya = floor(olddisp.at(a)/L);
				}
			else
				{
				xxa = olddisp.at(a)%L ;
				yya = floor(olddisp.at(a)/L) + 1;
				}

			// Posizione a-1
			if (a==0) {am1=N-1;}
			else {am1 = a-1;}

			if ((olddisp.at(am1))%L == 0)
				{
				xxam1 = L ;
				yyam1 = floor(olddisp.at(am1)/L);
				}
			else
				{
				xxam1 = olddisp.at(am1)%L ;
				yyam1 = floor(olddisp.at(am1)/L) + 1;
				}

			// Posizione a+1
			if (a==N-1) {ap1=0;}
			else {ap1 = a+1;}

			if ((olddisp.at(ap1))%L == 0)
				{
				xxap1 = L ;
				yyap1 = floor(olddisp.at(ap1)/L);
				}
			else
				{
				xxap1 = olddisp.at(ap1)%L ;
				yyap1 = floor(olddisp.at(ap1)/L) + 1;
				}

			// Posizione b
			if ((olddisp.at(b))%L == 0)
				{
				xxb = L ;
				yyb = floor(olddisp.at(b)/L);
				}
			else
				{
				xxb = olddisp.at(b)%L ;
				yyb = floor(olddisp.at(b)/L) + 1;
				}

			// Posizione b-1
			if (b==0) {bm1=N-1;}
			else {bm1 = b-1;}

			if ((olddisp.at(bm1))%L == 0)
				{
				xxbm1 = L ;
				yybm1 = floor(olddisp.at(bm1)/L);
				}
			else
				{
				xxbm1 = olddisp.at(bm1)%L ;
				yybm1 = floor(olddisp.at(bm1)/L) + 1;
				}

			// Posizione b+1
			if (b==N-1) {bp1=0;}
			else {bp1 = b+1;}

			if ((olddisp.at(bp1))%L == 0)
				{
				xxbp1 = L;
				yybp1 = floor(olddisp.at(bp1)/L);
				}
			else
				{
				xxbp1 = olddisp.at(bp1)%L ;
				yybp1 = floor(olddisp.at(bp1)/L) + 1;
				}

			// Nuova Confgurazione
			//______________________________________________________________
			if ((newdisp.at(a))%L == 0)
				{
				xa = L ;
				ya = floor(newdisp.at(a)/L);
				}
			else
				{
				xa = newdisp.at(a)%L ;
				ya = floor(newdisp.at(a)/L) + 1;
				}

			// Posizione a-1
			if (a==0) {am1=N-1;}
			else {am1 = a-1;}

			if ((newdisp.at(am1))%L == 0)
				{
				xam1 = L ;
				yam1 = floor(newdisp.at(am1)/L);
				}
			else
				{
				xam1 = newdisp.at(am1)%L ;
				yam1 = floor(newdisp.at(am1)/L) + 1;
				}

			// Posizione a+1
			if (a==N-1) {ap1=0;}
			else {ap1 = a+1;}

			if ((newdisp.at(ap1))%L == 0)
				{
				xap1 = L ;
				yap1 = floor(newdisp.at(ap1)/L);
				}
			else
				{
				xap1 = newdisp.at(ap1)%L ;
				yap1 = floor(newdisp.at(ap1)/L) + 1;
				}

			// Posizione b
			if ((newdisp.at(b))%L == 0)
				{
				xb = L ;
				yb = floor(newdisp.at(b)/L);
				}
			else
				{
				xb = newdisp.at(b)%L ;
				yb = floor(newdisp.at(b)/L) + 1;
				}

			// Posizione b-1
			if (b==0) {bm1=N-1;}
			else {bm1 = b-1;}

			if ((newdisp.at(bm1))%L == 0)
				{
				xbm1 = L ;
				ybm1 = floor(newdisp.at(bm1)/L);
				}
			else
				{
				xbm1 = newdisp.at(bm1)%L ;
				ybm1 = floor(newdisp.at(bm1)/L) + 1;
				}

			// Posizione b+1
			if (b==N-1) {bp1=0;}
			else {bp1 = b+1;}

			if ((newdisp.at(bp1))%L == 0)
				{
				xbp1 = L;
				ybp1 = floor(newdisp.at(bp1)/L);
				}
			else
				{
				xbp1 = newdisp.at(bp1)%L ;
				ybp1 = floor(newdisp.at(bp1)/L) + 1;
				}	
			//______________________________________________________________		

			// Calcolo il delta 
			Etrial = sqrt((xap1-xa)*(xap1-xa) + (yap1-ya)*(yap1-ya)) + sqrt((xam1-xa)*(xam1-xa) + (yam1-ya)*(yam1-ya)) + sqrt((xbp1-xb)*(xbp1-xb) + (ybp1-yb)*(ybp1-yb)) + sqrt((xbm1-xb)*(xbm1-xb) + (ybm1-yb)*(ybm1-yb));
			Eold = sqrt((xxap1-xxa)*(xxap1-xxa) + (yyap1-yya)*(yyap1-yya)) + sqrt((xxam1-xxa)*(xxam1-xxa) + (yyam1-yya)*(yyam1-yya)) + sqrt((xxbp1-xxb)*(xxbp1-xxb) + (yybp1-yyb)*(yybp1-yyb)) + sqrt((xxbm1-xxb)*(xxbm1-xxb) + (yybm1-yyb)*(yybm1-yyb));
			deltaE = Etrial -Eold;			

			// Testo di Metropolis
			if ( deltaE < 0 ) 
				{pacc = 1;}
			else 
				{pacc = 1/(exp(deltaE*Bl));}

			c = ran0(&idum); // Genero un numero casuale in [1,0]
			
			if (c < pacc ) // Se la probabilità di accettazione è maggiore del numero generato prendo la nuova configurazione e la sua energia
				{
				olddisp = newdisp;	
				Ej = Ej + deltaE ;
				Ej2 = Ej2 + (deltaE)*(deltaE);
				cont = cont +1;	// Conto le energie accettate ogni simulazione		
				}
			}

		for (int i = 0; i < N ; i++) // Calcolo le posizioni delle città nella griglia per l'ultima configurazione accettata
				{
				if ((olddisp.at(i))%L == 0)
					{
					x.at(i) = L ;
					y.at(i) = floor(olddisp.at(i)/L);
					}
				else
					{
					x.at(i) = olddisp.at(i)%L ;
					y.at(i) = floor(olddisp.at(i)/L) + 1;
					}
				}

		Ek1=0;
			
		for (int i = 0; i < N-1; i++) // Calcolo l'energia della nuova configurazione come somma delle distanze tra città successive
			{
			Ek1 += sqrt((x.at(i)-x.at(i+1))*(x.at(i)-x.at(i+1))+(y.at(i)-y.at(i+1))*(y.at(i)-y.at(i+1)));
			}
				
		Ek1 += sqrt((x.at(N-1)-x.at(0))*(x.at(N-1)-x.at(0))+(y.at(N-1)-y.at(0))*(y.at(N-1)-y.at(0))); // Aggiungo la distanza per tornare dall'ultima città alla prima
	
		ak = cont*1.0L/(N*Nk); // Calcolo il tasso di accettazione per ogni set di simulazioni

		grafen << passo << " " << Ek1 << endl;
		cout << "Passo: " << passo << endl << "Varianza: " << varEk << endl << "Beta: " << Bl << endl << "Tasso di accettazione: " << ak << endl<< "Energia accettata: " << Ek1 << endl;

		passo++;	
		cout << "-------------------------------------------" << endl;	
		}

	cout << "Energia minima: " << Ek1 << endl; 

	cout << "Configurazione iniziale: ";
	for (int i = 0; i < N-1 ; i++)
	{cout << confiniziale.at(i) << ", ";}
	cout << confiniziale.at(N-1) << endl;

	cout << "Configurazione finale: ";
	for (int i = 0; i < N-1 ; i++)
	{cout << olddisp.at(i) << ", ";}
	cout << olddisp.at(N-1) << endl;

	for (int i = 0; i < N ; i++) // Calcolo le posizioni delle città nella griglia
		{
		if ((olddisp.at(i))%L == 0)
			{
			x.at(i) = L ;
			y.at(i) = floor(olddisp.at(i)/L);
			}
		else
			{
			x.at(i) = olddisp.at(i)%L ;
			y.at(i) = floor(olddisp.at(i)/L) + 1;
			}
		}

	for (int i = 0; i < N ; i++)
	{periniz << xstart.at(i) << " " << ystart.at(i) << endl;}

	for (int i = 0; i < N ; i++)
	{perfin << x.at(i) << " " << y.at(i) << endl;}

	return 0;
}

double ran0(long *idum) // Generatore di numeri casuali di tipo LCG con distribuzine uniforme in [0,1]
{
	double ans;

  *idum=(IA*(*idum))%IM;  // Aggiornamento di idum
  ans=AM*(*idum); // Conversione di idum in un numero uniforme in [0,1]

  return ans;
}

