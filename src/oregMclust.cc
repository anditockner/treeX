/*
 * Circular/linear regression clustering via redescending M-estimators.
 *
 * Vendored from the CRAN package 'edci' version 1.1-3 (GPL-2)
 * by Tim Garlipp, TimGarlipp@gmx.de (see https://CRAN.R-project.org/package=edci).
 * Previously shipped in treeX as a precompiled binary ('edc'); included
 * here as source so the package compiles natively on all platforms
 * (including macOS on Intel and Apple Silicon).
 */
#include <cmath>
#include <cassert>
#include <R.h>

/* R >= 4.5 enables STRICT_R_HEADERS by default, so the legacy PI macro from
   R_ext/Constants.h is no longer provided - define it ourselves. */
#ifndef PI
#define PI 3.14159265358979323846
#endif
using namespace std;

#define QUAD(x) ((x)*(x))

#define BFAK 1.6 // Faktor zum Beschleunigen der Schrittweite
#define DFAK 0.7 // Faktor zum D?mpfen der Schrittweite (Armijo-Schrittweite)
#define AFAK 0.5 // Faktor f?r die Abbruchbedingun im Armijo-Verfahren

//#define PI 3.141592654

// Negative Dichte der Std.Normalverteilung und Ableitungen
double rho (double x){
    return -exp(-(x*x)/2)/2.506628;
}

double rho1 (double x){
    return -rho(x)*x;
}

double rho2 (double x){
    return -(x*rho1(x)+rho(x));
}

// Vektor der obigen Werte (schneller)
void rhov (double x, double *r){
    r[0] = -exp(-(x*x)/2)/2.506628;
    r[1] = -r[0]*x;
    r[2] = -(x*r[1]+r[0]);

    return;
}


// Orthogonale Regression
extern "C" {
    
    void c_oregMclust (
	double *const x,           // x- und y- Koordinaten der Beobachtungen
	double *const y, 
	int    *const n,           // Anzahl der Beobachtungen
	double *const bandbreite,
	int    *const method,      // Methode der Bestimmung der Startwerte:
	                           // 0: F?r jede Beobachtung 'anzwinkel'
 	                           //    verschiedene Winkel
        	                   // 1: F?r jede Beobachtung den Winkel
       	                           //    aus 'startwinkel'
	                           // 2: Nur die eine Gerade mit Parametern
	                           //    'salpha', 'sbeta'
	                           // 3: Alle geraden durch jeweils zwei
	                           //    Beobachtungen
	                           // Die jeweils nicht ben?tigten Parameter
	                           // d?rfen 'NULL' sein.
	int    *const anzwinkel,   // Anzahl der Startwerte f?r den Winkel
	double *const startwinkel, // Startwerte f?r die Winkel
	double *const salpha,      // Parameter der Startgeraden
	double *const sbeta,
	int    *const brmaxit,     // Abbruch nach brmaxit Iterationen
	// R?ckgabewerte
	double *const alpha,       // Parameter der gefundenen Geraden
	double *const beta, 
	double *const value,       // und zugehoehriger Wert der scorefkt.
	int          *count)       // Anzahl der gefundenen Geraden
	//double *const xp,          // orth. Projektion der Beobachtungen auf 
	//double *const yp)          // die gefundenen Geraden bei Methode 1
	                           // oder Methode 0 mit anzwinkel=1
	                           
    {
	int i, j, k;
	double a, b;                       // Schaetzungen f?r alpha und beta
	double diff_a, diff_b;             // Schrittl?ngen
	double h, ha, hb, haa, hab, hbb;   // Wert von Hn sowie erste und 
	                                   // zweite Ableitungen
	double det;                        // Determinante der Hessematrix
	double hilf, hilfa, hilfb, hilfaa; // Innerer Term + innere Ableitungen
	double hilf_h;                     // Wert von Hn an anderer Stelle
	double fak;                        // Faktor zum Beschleunigen oder
	                                   // D?mpfen der Schrittweite
	double r[3];                       // Vektor der Scorefkt. und deren
	                                   // Ableitungen
	double sn = *bandbreite;           // Bandbreite
	int anz;
	int begrenzen;
	double altdiff_a;
	double altdiff_b;

	int anzstartpunkte;
	int anzstartwinkel;

	assert (0 <= (*method) && (*method) <= 3);

	if (*method == 2){
	    anzstartpunkte = 1;
	    anzstartwinkel = 1;
	}
	else {
	    anzstartpunkte = *n;
	    if (*method == 0)
		anzstartwinkel = *anzwinkel;
	    else if (*method == 1)
		anzstartwinkel = 1;
	    else if (*method == 3)
		anzstartwinkel = *n;
	    else 
		anzstartwinkel = 0; // Nur, um den Compiler zu beruhigen
	}

	*count = 0;
        
	for (i=0; i < anzstartpunkte; i++){
	    for (j=0; j < anzstartwinkel; j++){
		if ((*method!=2) && (*count % 100 == 0))
		    Rprintf("starting point: %i/%i\n",
			    *count,anzstartpunkte*anzstartwinkel);
		if ((*method!=3) || 
		    ((j>i) && ( (x[j] != x[i]) || (y[j] != y[i]))))
		{
		    if (*method==2){
			a = *salpha;
			b = *sbeta;
		    }
		    else {
			if (*method==0)
			    a = (double) j/(*anzwinkel)*PI;
			else if (*method==1)
			    a = startwinkel[i];
			else if (*method==3)
			    // Winkel der Gerade durch die Punkte (x[i],y[i]) 
			    // und  (x[j],y[j]) bestimmen
			    if (y[j]==y[i])
				a = PI/2;
			    else
				a = atan(-(x[j]-x[i])/(y[j]-y[i]));
			else
			    a = 0; // Nur, um den Compiler zu beruhigen
			b = cos(a)*x[i]+sin(a)*y[i];
		    }
		
  		    anz = 0;
		    begrenzen = 1;
		    altdiff_a = 0;
		    altdiff_b = 0;
		    do {
			//Werte und Ableitungen von Hn bestimmen
			h   = 0; 
			ha  = hb  = 0; 
			haa = hab = hbb = 0;
			for (k=0; k < *n; k++){
			    // Inneren Term + innere Ableitungen berechnen
			    hilf   = 1/sn*(cos(a)*x[k]+sin(a)*y[k]-b);
			    hilfa  = 1/sn*(-sin(a)*x[k]+cos(a)*y[k]);
			    hilfb  = -1/sn;
			    hilfaa = 1/sn*(-cos(a)*x[k]-sin(a)*y[k]);
			    
			    rhov (hilf, r);

			    h   = h   + r[0];
			    ha  = ha  + r[1]*hilfa;
			    hb  = hb  + r[1]*hilfb;
			    haa = haa + r[2]*hilfa*hilfa+r[1]*hilfaa;
			    hab = hab + r[2]*hilfb*hilfa;
			    hbb = hbb + r[2]*hilfb*hilfb;

			}
			
			det = haa*hbb - hab*hab;
			
			if (det > 0 && haa > 0){             // Newtown
			    // - (Inverse der Hessematrix) * Gradient
			    diff_a = -(hbb/det*ha-hab/det*hb);
			    diff_b = -(-hab/det*ha+haa/det*hb);

			    fak = DFAK; // Der Newtonschritt mu? noch ged?mpft
			                // werden
      			}
			else {                               // Stepest Descent
			    // Die Richtung ist der Gradient
			    diff_a = -ha;
			    diff_b = -hb;

			    // Beschleunigung des Verfahrens durch
			    // Vergr??erung der Schrittweite:
			    // Solange Hn in Richtung des Gradienten kleiner 
			    // wird, kann die Schrittweite vergr??ert werden

			    fak = 1/BFAK;
			    do {
				fak = fak * BFAK;
				// Hn(a+diff_a, b+diff_b) bestimmen
				hilf_h = 0;
				for (k=0; k < *n; k++)
				    hilf_h = hilf_h + 
					rho(1/sn*(cos(a+fak*diff_a)*x[k]+
						  sin(a+fak*diff_a)*y[k] - 
						  (b+fak*diff_b)));
			    } while (h - hilf_h > 0.00001);
			    fak = fak / BFAK; // Sollte die Schleife schon 
			                      // nach einem Druchlauf
			                      // abgebrochen werden, ist fak 
			                      // jetzt < 1 und die
			                      // Schrittweite wird gleich 
			                      // ged?mpft.
			}

			// D?mpfung der Schrittweite nach 'Armijo's Rule', 
			// falls nicht gerade beschleunigt wurde

			if (fak < 1){
			    fak  = fak/DFAK;
			    do {
				fak = fak * DFAK;
				// Hn(a+diff_a, b+diff_b) bestimmen
				hilf_h = 0;
				for (k=0; k < *n; k++)
				    hilf_h = hilf_h + 
					rho(1/sn*(cos(a+fak*diff_a)*x[k]+
						  sin(a+fak*diff_a)*y[k] - 
						  (b+fak*diff_b)));
			    } while ((hilf_h - (h + AFAK*fak * 
						(ha*diff_a+hb*diff_b)) 
				      > 0.00001) && 
				     ((begrenzen == 1) | (fak > 0.0001)));
			}

			if ((altdiff_a*diff_a < 0) | (altdiff_b*diff_b < 0))
			    begrenzen = 0;
			altdiff_a = diff_a;
			altdiff_b = diff_b;

			a = a + fak*diff_a;
			b = b + fak*diff_b;

			anz ++;
			
		    } while (fak*(fabs(diff_a)+fabs(diff_b)) > 0.00001 && 
			     anz <= (*brmaxit));

		    if (anz <= (*brmaxit)) {
			alpha[*count] = a;
			beta[*count]  = b;
			value[*count] = h;
			*count = *count + 1;
		    }
		}
	    }
	}
	
	return;
    } // c_oregMclust

    //*****************************************************************************************************

    // Orthogonale Regression f?r Kreise
    void c_oregMcirc (double *const x, double *const y,  // x- und y- Koordinaten der Beobachtungen
		      int    *const n,                   // Anzahl der Beobachtungen
		      double *const bandbreite,
		      int    *const method,              // 0: alle, 1: einer mit Kreis durch sij
		                                         // 2: scount mit Param. scx, scy, sr 
		      int    *const si1,                 // Startpunkte
		      int    *const si2,
		      int    *const si3,
		      double *const maxstartrad,
		      double *const startcx,
		      double *const startcy,
		      double *const startr,
		      int    *const startcount,
		      double *const breakminx,           // Wenn die Parameter diese Grenzen
		      double *const breakmaxx,           // verlassen oder die Anzahl der
		      double *const breakminy,           // Minimierungsiterationen "breakmaxit"
		      double *const breakmaxy,           // ?bersteigt, wird abgebrochen
		      double *const breakminr,
		      double *const breakmaxr,
		      int    *const breakmaxit,
		      // R?ckgabewerte
		      double *const cx,                  // Mittelpunkte und
		      double *const cy,
		      double *const r,                   // Radien der gefundenen Kreise
		      double *const value,               // und zugeh?rige Werte der Scorefunktion
		      int          *count)               // Anzahl der gefundenen Kreise
    {
	int i1, i2, i3, k;

	double a1, a2, b;                   // Schaetzungen f?r cx/cy und r
	double diff_a1, diff_a2, diff_b;    // Schrittl?ngen
	double h, ha1, ha2, hb,             // Wert von Hn sowie erste und zweite Ableitungen
	    ha1a1, ha1a2, ha1b,
	    ha2a1, ha2a2, ha2b,
	    hba1,  hba2,  hbb;  
	double det, det2;                   // Determinante und 2x2 Unterdeterminante der Hessematrix
	double in, hilf, hilfa1, hilfa2,    // Innerer Term + innere Ableitungen
	    hilfb, hilfa1a1, hilfa2a2, hilfaa;
	double hilf_h, merk_h ;             // Wert von Hn an anderer Stelle
	double fak;                         // Faktor zum Beschleunigen oder D?mpfen der Schrittweite
	double rhovec[3];                   // Vektor der Scorefkt. und deren Ableitungen
	double sn = *bandbreite;            // Bandbreite
	int anz;
	int begrenzen;
	double altdiff_a1, altdiff_a2;
	double altdiff_b;

	double x1,x2,x3, y1,y2,y3;          // Punkte zur Berechnung der Startkreise
	double mb, mc, bb, bc;              // Steigungen und Achsenabschnitte der Mittelsenkrechten

	anz = (*n)*((*n)-1)/2;
	int    abbr;
	double breakmindiff, sbmd;          // Minimierung abbrechen, wenn Parameter sich um weniger
	                                    // als breakmindiff ver?ndern

	double startx, starty, startr1;

	int    von,bis;

	int    br;

	*count = 0;
	breakmindiff = 0.00001*(((*breakmaxx)-(*breakminx))>((*breakmaxy)-(*breakminy))?
				((*breakmaxx)-(*breakminx)):((*breakmaxy)-(*breakminy)));
	sbmd = 0.0001*sqrt(((*breakmaxx)-(*breakminx))>((*breakmaxy)-(*breakminy))?
			   ((*breakmaxx)-(*breakminx)):((*breakmaxy)-(*breakminy)));

	// F?r je drei Punkte den durch diese laufenden Kreis als Startwert

	von = 0;
	if (*method==0)
	    bis = *n;
	else if (*method==1)
	    bis = 1;
	else if (*method==2)
	    bis = *startcount;
	else
	    return;
	
	for (i1=von; i1 < bis; i1++){
	    if (*method==0)
		von=i1+1;
	    else if (*method==2)
		von=bis-1;
	    for (i2=von; i2 < bis; i2++) {
		if (*method==0)
		    von=i2+1;
		for (i3=von; i3 < bis; i3++) {
		    if (*method==2){
			a1 = startcx[i1];
			a2 = startcy[i1];
			b  = startr[i1];
			abbr = 0;
		    }
		    else {
			if (*method==1){
			    i1 = *si1;
			    i2 = *si2;
			    i3 = *si3;
			}
			if (*method==0)
			    Rprintf("starting points %i/%i/%i of %i:",i1,i2,i3,*n);
			
			// Die Punkte so permutieren, da? y3!=y1!=y2, falls die drei Punkte 
			// nicht auf einer geraden liegen
			if (y[i1] == y[i2]){
			    x1 = x[i3]; y1 = y[i3];
			    x2 = x[i2]; y2 = y[i2];
			    x3 = x[i1]; y3 = y[i1];
			}
			else if (y[i1] == y[i3]){
			    x1 = x[i2]; y1 = y[i2];
			    x2 = x[i1]; y2 = y[i1];
			    x3 = x[i3]; y3 = y[i3];
			}
			else {
			    x1 = x[i1]; y1 = y[i1];
			    x2 = x[i2]; y2 = y[i2];
			    x3 = x[i3]; y3 = y[i3];
			}
			
			// Nur wenn die drei Punkte nicht auf einer Geraden liegen, bilden sie
			// einen Startwert
			if ((y1 != y2 || y1 != y3) &&
			    (x1 != x2 || x1 != x3 || x2 != x3)){  

			    abbr = 0;

			    mc = - (x2-x1)/(y2-y1);         // Steigung und Achsenabschnitt der 
			                                    // Mittelsenkrechten zwischen
			    bc = (y2+y1)/2 - mc*(x2+x1)/2;  // (x1,y1) und (x2,y2)
			
			    mb = - (x3-x1)/(y3-y1);         // Dito f?r (x1,y1) und (x3,y3)
			    bb = (y3+y1)/2 - mb*(x3+x1)/2;
			
			    a1 = (bb-bc)/(mc-mb);       // Schnittpunkt der beiden Mittelsenkrechten
			    a2 = mc*a1+bc;              // ist der Mittelpunkt des Umkreises
			
			    b = sqrt(QUAD(a1-x1)+QUAD(a2-y1));  // Radius als Abstand eines Punktes
			                                        // zum Mittelpunkt
			}
			else {
			    abbr = 1;
			    a1   = 0; // Compiler beruhigen
			    a2   = 0;
			    b    = 0;
			}
		    }

		    // Falls nicht, eigentliche Minimierung starten
		    if (abbr == 0 && b<=(*maxstartrad)){
		        anz        = 0;
			begrenzen  = 1;
			altdiff_a1 = 0;
			altdiff_a2 = 0;
			altdiff_b  = 0;
			startx     = a1;
			starty     = a2;
			startr1    = b;

			br = 0;
			do {
			    //Werte und Ableitungen von Hn bestimmen
			    h     = 0; 
			    ha1   = ha2   = hb   = 0; 
			    ha1a1 = ha1a2 = ha1b = 0;
			    ha2a1 = ha2a2 = ha2b = 0;
			    hba1  = hba2  = hbb  = 0;
			    for (k=0; k < *n; k++){
				// Inneren Term + innere Ableitungen berechnen
				in       = QUAD(x[k]-a1)+QUAD(y[k]-a2);
				hilf     = (double) 1/sn*(sqrt(in)-b);
				hilfa1   = (double) -1/sn*(x[k]-a1)*1/sqrt(in);
				hilfa2   = (double) -1/sn*(y[k]-a2)*1/sqrt(in);
				hilfb    = (double) -1/sn;
				hilfa1a1 = (double) -1/sn*(1/(in*sqrt(in))*QUAD(x[k]-a1)-1/sqrt(in));
				hilfa2a2 = (double) -1/sn*(1/(in*sqrt(in))*QUAD(y[k]-a2)-1/sqrt(in));
				hilfaa   = (double) -1/sn*1/(in*sqrt(in))*(x[k]-a1)*(y[k]-a2);
				    
				rhov (hilf, rhovec);
				
				h   = h   + rhovec[0];
				    
				// erste Ableitungen
				ha1 = ha1 + rhovec[1]*hilfa1;
				ha2 = ha2 + rhovec[1]*hilfa2;
				hb  = hb  + rhovec[1]*hilfb;
				
				// zweite Ableitungen
				ha1a1 = ha1a1 + rhovec[2]*hilfa1 * hilfa1  +  rhovec[1]*hilfa1a1;
				ha1a2 = ha1a2 + rhovec[2]*hilfa2 * hilfa1  +  rhovec[1]*hilfaa;
				ha1b  = ha1b  + rhovec[2]*hilfb  * hilfa2;
				
				ha2a1 = ha2a1 + rhovec[2]*hilfa1 * hilfa2  +  rhovec[1]*hilfaa;
				ha2a2 = ha2a2 + rhovec[2]*hilfa2 * hilfa2  +  rhovec[1]*hilfa2a2;
				ha2b  = ha2b  + rhovec[2]*hilfb  * hilfa2;
				
				hba1  = hba1  + rhovec[2]*hilfa1 * hilfb;
				hba2  = hba2  + rhovec[2]*hilfa2 * hilfb;
				hbb   = hbb   + rhovec[2]*hilfb  * hilfb;
			    }
			    
			    det = ha1a1*ha2a2*hbb + ha1a2*ha2b*hba1 + ha1b*ha2a1*hba2 -
				( hba1*ha2a2*ha1b + hba2*ha2b*ha1a1 + hbb*ha2a1*ha1a2);
			    
			    det2 = ha1a1*ha2a2 - ha1a2*ha2a1;
			    
			    if (det > 0 && det2 > 0 && ha1a1 > 0){            // Newtown
				// - (Inverse der Hessematrix) * Gradient
				diff_a1 = -(( -ha2b*hba2  +ha2a2*hbb)  /det*ha1 
					    +( ha1b*hba2  -ha2a2*hbb)  /det*ha2 
					    +(-ha1b*ha2a2 +ha1a1*hbb)  /det*hb);
				diff_a2 = -((  ha2a1*hba1 -ha2a1*hbb)  /det*ha1 
					    +(-ha1b*hba1  +ha1a1*hbb)  /det*ha2 
					    +( ha1b*ha2a1 -ha1a1*ha2b) /det*hb);
				diff_b  = -(( -ha2a2*hba1 +ha2a1*hba2) /det*ha1 
						+( ha1a1*hba1 -ha1a1*hba2) /det*ha2 
						+(-ha2a2*ha2a1+ha1a1*ha2a2)/det*hb);
				    
				fak = DFAK; // Der Newtonschritt mu? noch ged?mpft werden
			    }
			    else {                                           // Stepest Descent
				// Die Richtung ist der Gradient
				diff_a1 = -ha1;
				diff_a2 = -ha2;
				diff_b  = -hb;
				
				// Beschleunigung des Verfahrens durch
				// Vergr??erung der Schrittweite:
				// Solange Hn in Richtung des Gradienten kleiner wird,
				// kann die Schrittweite vergr??ert werden
				
				fak = 1;
				if (begrenzen == 1){

				    fak    = 1/BFAK;
				    hilf_h = h;
				    do {
					merk_h=hilf_h;
					fak = fak * BFAK;
					// Hn(a1+diff_a1, a2+diff_a2, b+diff_b) bestimmen
					hilf_h = 0;
					for (k=0; k < *n; k++)
					    hilf_h = hilf_h + 
						rho(1/sn*(sqrt(QUAD(x[k]-(a1+fak*diff_a1))+
							       QUAD(y[k]-(a2+fak*diff_a2)))-
							  (b+fak*diff_b)));
				    } while (merk_h - hilf_h > 0.1*sbmd);
				}
				fak = fak / BFAK; // Sollte die Schleife schon nach einem 
				                  // Durchlauf abgebrochen werden, ist fak 
				                  // jetzt < 1 und die Schrittweite 
				                  // wird gleich ged?mpft.
			    }
			    
			    // D?mpfung der Schrittweite nach 'Armijo's Rule', falls nicht gerade
			    // Beschleunigt wurde
			    
			    if (fak < 1){
				fak  = fak/DFAK;
				do {
				    fak = fak * DFAK;
				    // Hn(a1+diff_a1, a2+diff_a2, b+diff_b) bestimmen
				    hilf_h = 0;
				    for (k=0; k < *n; k++)
					hilf_h = hilf_h + 
					    rho(1/sn*(sqrt(QUAD(x[k]-(a1+fak*diff_a1))+
							   QUAD(y[k]-(a2+fak*diff_a2)))-
						      (b+fak*diff_b)));
				} while ((hilf_h - (h + AFAK*fak * 
						    (ha1*diff_a1+ha2*diff_a2+hb*diff_b)) 
					  > 0.01*sbmd)
					 && ((begrenzen == 1) | (fak > 0.1*sbmd)));
			    }
			    
			    if ((altdiff_a1*diff_a1 < 0) | (altdiff_a2*diff_a2 < 0) |
				(altdiff_b*diff_b < 0)){
				begrenzen = 0;
			    }
			    altdiff_a1 = diff_a1;
			    altdiff_a2 = diff_a2;
			    altdiff_b  = diff_b;
			    
			    a1 = a1 + fak*diff_a1;
			    a2 = a2 + fak*diff_a2;
			    b  = b  + fak*diff_b;
				
			    anz ++;

			    if ( (*breakminx) > a1 || a1 > (*breakmaxx) ||
				 (*breakminy) > a2 || a2 > (*breakmaxy) ||
				 (*breakminr) > b  || b  > (*breakmaxr) ||
				 (anz > (*breakmaxit)))
				br=1;
			    
			} while (fak*(fabs(diff_a1)+fabs(diff_a2)+fabs(diff_b)) >= breakmindiff
				 && br==0);

 			if (br==0 && a1+a2+b+h==a1+a2+b+h && 
			    (int)(h*1000000)!=0)
			{
			    cx[*count]    = a1;
			    cy[*count]    = a2;
			    r[*count]     = b;
			    value[*count] = h;
			    *count = *count + 1;
			}
		    }
		}
	    }
	}
	
	return;
    }
    

} // extern
