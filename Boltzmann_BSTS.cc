/*
 *  BOLTZMANN EQUATION SOLVED FOR PHONON MEDIATED EXCITED STATE DYNAMICS IN BSTS TOPOLOGICAL INSULATORS
 *  2018, O. Abdurazakov, oabdura@ncsu.edu 
 * 
 * This code solves the Boltzmann Transport Equation in the presence of electron-phonon scattering for a linear/square electronic density of states relevant to BSTS topological insulators.
 * Output: Electron/phonon density distributions as a function of time: f(t)/n(t)   
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include <assert.h>
#include <iomanip>
    using namespace std;

#include "matrix.h"

#define band_bottom 0.0		// Bottom of the conduction band 
#define band_top    2.0		// Top of the conduction band 
#define phonon_band_top 0.02	// Max phonon energy
#define dw 0.001		// Energy step 

#define lambda 70.00		// Electron-lattice coupling constant
#define G_escape (0.000)	//Escape rate from the conduction band 

#define tmin 0.0		// Min time
#define tmax 2000.00		// Max time
#define dt   0.01		// Time step

#define kb (8.617e-5)		// 1 Kelvin = 8.617e-5 eV
#define phonon_T   (100.0 * kb) // Lattice temperature 
#define electron_T 0.2		// Temperature of the initial hot (Boltzmann) distribution
#define f_init 0.2		// Prefactor of the initial hot distribution
#define occ 1 			// Phase space constraint (i.e., Pauli blocking) imposed if 1

#define update_phonons false	// Lattice temperature stays fixed if false 

int Nw = int((band_top-band_bottom)/dw);
int Nw_ph = int((phonon_band_top)/dw);
int Nt = int((tmax-tmin)/dt);

double fermi(double x,double T)
{
    return 1./(exp(x/T)+1.0);
}

double bose(double x,double T)
{
    return 1./(exp(x/T)-1.0);
}

void setup_electron_dos(vector<double> &electron_dos)
{
    // Set up electron dos
    ofstream edfile("electron_dos.dat");
    electron_dos[0] = 0.0;
    for(int iw=1; iw < Nw; iw++)
    {
        double w = band_bottom + iw * dw;
        //electron_dos[iw] = 1.0 + (w-band_top/2)/2;
        //electron_dos[iw] = sqrt(w);
        electron_dos[iw] = w;
        edfile << iw << " " << w << " " << electron_dos[iw] << endl;
    }
    edfile.close();
}

void setup_phonon_dos_gauss(vector<double> &phonon_dos)
{
    double mu = 0.010;
    double sigma = 0.001;

    // Set up phonon dos
    ofstream pdfile("phonon_dos.dat");
    for(int iw=1; iw < Nw_ph; iw++)
    {
        double w = iw * dw;
        phonon_dos[iw] = exp(-0.5*(w-mu)*(w-mu)/(sigma*sigma)) ;
        pdfile << iw << " " << w << " " << phonon_dos[iw] << endl;
    }
}

void setup_phonon_dos(vector<double> &phonon_dos)
{
    // Set up phonon dos
    ofstream pdfile("phonon_dos.dat");
    for(int iw=0; iw < Nw_ph; iw++)
    {
        double w = iw * dw ;

        if(w < 0.004)
            phonon_dos[iw] = 6.25*1e4*w*w;
        else if(w >=0.004 && w < 0.018)
            phonon_dos[iw]=1.0;
        else
            phonon_dos[iw] = 0.0;

        pdfile << iw << " " << w << " " << phonon_dos[iw] << endl;
    }
}

void dump_f(matrix<double> &f)
{
    ofstream outfile("f.dat");
    for(int it=0; it < Nt; it++)
    {
        for(int iw=0; iw < Nw; iw++)
        {
            	outfile << it << " " << iw << " " << it*dt << " " << iw*dw << " " << f(iw,it) << endl;
        }
        outfile << endl;
    }
    outfile.close();
}


void collision(double* f, double* n, double* coll, vector<double> &electron_dos, vector<double> &phonon_dos)
{
    #pragma omp parallel for
    for(int ix=0; ix < Nw; ix++) // ix is electron
    {
        double x = band_bottom + dw*ix;


        for(int iw=1; iw < Nw_ph; iw++) // iw is phonon
        {
            double fa=0, fe=0, ee=0, ea=0;
            double w = iw * dw ;

            if(phonon_dos[iw] < 1e-10)
            {
                coll[ix + Nw*iw] = 0.;
                continue;
            }

            if(x-w > band_bottom)
            {
                fa +=      n[iw]  * f[ix-iw] * (1. - occ * f[ix]   ) * electron_dos[ix-iw] * phonon_dos[iw];
                ee += (1 + n[iw]) * f[ix]    * (1. - occ * f[ix-iw]) * electron_dos[ix-iw] * phonon_dos[iw];
                assert(ix-iw >= 0);
            }
            if(x+w < band_top)
            {
                fe += (1 + n[iw]) * f[ix+iw] * (1. - occ * f[ix]   ) * electron_dos[ix+iw] * phonon_dos[iw];
                ea +=      n[iw]  * f[ix]    * (1. - occ * f[ix+iw]) * electron_dos[ix+iw] * phonon_dos[iw];
                assert(ix+iw < Nw);
            }

            coll[ix + Nw*iw] = lambda * (fa + fe - ea - ee) * dw ;
        }
    }
}

int main()
{
    vector<double> electron_dos(Nw);
    vector<double> phonon_dos(Nw_ph);

    setup_electron_dos(electron_dos);
    setup_phonon_dos(phonon_dos);
    //setup_phonon_dos_gauss(phonon_dos);

    matrix<double> f(Nw,Nt);
    matrix<double> n(Nw_ph,Nt);
    double* coll = new double[Nw*Nw_ph];
    double* collnext = new double[Nw*Nw_ph];

    ofstream outfile("f.dat");
    ofstream out2file("n.dat");

    // Setup initial distribution
    {
        int it=0;
        for(int iw=0; iw < Nw; iw++)
        {
            double w = band_bottom + iw * dw;
            //f(iw,it) = fermi(w, phonon_T);
            f(iw,it) = f_init * exp(-w/electron_T);
	    if((it+1)%100==0 && iw*dw + band_bottom < 1)
            	outfile << it << " " << iw << " " << it*dt << " " << iw*dw << " " << f(iw,it) << " " << 0. << " " << 0. << endl;
        }
        outfile << endl;

        for(int iwph=1; iwph < Nw_ph; iwph++)
        {
            double w = iwph * dw;
            n(iwph,it) = bose(w, phonon_T);
        }
    }

    // Calculate the it=0 collision integrals
    collision(f.col(0), n.col(0), coll, electron_dos, phonon_dos);

    double eloss =0 ;
    double toten, totenprev, toten0;
    double eph0, eph, ephprev;


    toten0 = 0;
    for(int iw=0; iw < Nw; iw++)
        toten0 += (band_bottom + iw * dw) * f(iw,0) * dw * electron_dos[iw];
    toten = toten0;
    eph0 = 0;
    for(int iwph=1; iwph < Nw_ph; iwph++)
        eph0 += n(iwph, 0) * dw * (iwph * dw) * phonon_dos[iwph];
    eph = eph0;

    /// MAIN TIME LOOP ===========================================================
    for(int it=0; it < Nt-1; it++)
    {

        // Expand n
        if(it==0)
            for(int iwph = 1; iwph < Nw_ph; iwph++)
                n(iwph, it+1) = n(iwph, it) ;
        else
            for(int iwph=1; iwph < Nw_ph; iwph++)
                n(iwph, it+1) = 2*n(iwph, it) - n(iwph, it-1);

        // Step based on old collision integral
        for(int ix=0; ix < Nw; ix++)
        {
            f(ix,it+1) = f(ix,it);
            for(int iwph = 1; iwph < Nw_ph; iwph++)
                f(ix,it+1) += dt * coll[ix + Nw * iwph];

            f(ix,it+1) -= dt * G_escape * f(ix,it);
        }


        // Compute the collision integral of the next step
        collision(f.col(it+1), n.col(it+1), collnext, electron_dos, phonon_dos);

        double fsum=0, nsum=0;
        totenprev = toten;
        toten = 0;

        // Step electrons
        for(int iw=0; iw < Nw; iw++)
        {
            // Step based on average of two collision integrals
            f(iw,it+1) = f(iw,it);
            for(int iwph = 1; iwph < Nw_ph; iwph++)
                f(iw,it+1) += 0.5 * dt * (coll[iw + Nw * iwph] + collnext[iw + Nw * iwph]);

            f(iw,it+1) -= dt * G_escape * f(iw,it);


            toten += (band_bottom + iw * dw) * f(iw,it+1) * dw * electron_dos[iw];





            // Sum tracker
            fsum += f(iw,it+1) * dw * electron_dos[iw];

            // Some output
            {
            if((it+1)%100==0 && iw*dw + band_bottom < 1)
                outfile << (it+1) << " " << iw << " " << (it+1)*dt << " " << iw*dw << " " << f(iw,it+1) << " " << f(iw,it+1)-f(iw,0) << " " << coll[iw] << endl;
            // Some checks
            if(f(iw,it+1) < 0)
            {
                //outfile.flush();
                //cout << iw << " " << it+1 << endl;
                //assert(false && "Zero at " && it+1 && " " && iw);
                if(band_bottom + iw * dw < 0.5)
                    cout <<  "Negative at time " << (it+1)*dt << " and energy " << band_bottom + iw*dw << endl;
                f(iw,it+1) = 0.;
            }
            }
        }


        // Compute energy lost per phonon mode and total energy in the phonons
        ephprev = eph;
        eph = 0.;
        double eloss=0;
        nsum = 0;
        for(int iwph = 1; iwph < Nw_ph; iwph++)
        {
            //cout << iwph << " " << Nw_ph << endl;
            double w = iwph * dw;
            double csum=0;
            for(int ix=0; ix < Nw; ix++)
            {
                double x = band_bottom + dw * ix;
                csum += dt * x * 0.5 * (coll[ix + Nw * iwph] + collnext[ix + Nw * iwph]) * electron_dos[ix];

            }
            //cout << "Lost to w= " << w << " : " << csum << endl;
            eloss += csum * dw;

            if(phonon_dos[iwph] > 1e-10)
            {
                //cout << "Updating w= " << iwph * dw << " " << dw * w * phonon_dos[iwph] * csum  / (dw * w * phonon_dos[iwph]) << endl;
                if(update_phonons)
                    n(iwph, it+1) = n(iwph, it) - csum / (w * phonon_dos[iwph]);
                else
                    n(iwph, it+1) = n(iwph, it);

                nsum += n(iwph, it+1);
            }

            if((it+1)%100==0)
            {
                out2file << (it+1) << " " << iwph << " " << (it+1)*dt << " " << iwph*dw << " " << n(iwph,it+1) << " " << n(iwph,it+1)-n(iwph,0) << endl;
            }

            eph += n(iwph, it+1) * dw * w * phonon_dos[iwph];
        }


        // Swap
        swap(coll, collnext);

       
        if((it+1)%100==0)
        {
            cout << setw(5) << setprecision(5) << (it+1)*dt << " " << setw(10) << setprecision(8) 
                << fsum << " " << nsum << " "
                /*
                << "\t|\t" << toten << " " << toten-totenprev << " " << toten-toten0 
                << "\t|\t" << eph << " " << eph-ephprev << " " << eph - eph0 
                << "\t|\t" << eloss
                */
                << "\t|\t" << toten << " " << eph 
                << "\t|\t" << toten-totenprev << " " << eph-ephprev 
                << "\t|\t" << toten-toten0 << " " << eph-eph0
                << "\t|\t" << eloss 
                << endl;
            outfile << endl;
            out2file << endl;
        }
    }

    delete[] coll;
    delete[] collnext;

    //dump_f(f);
}
