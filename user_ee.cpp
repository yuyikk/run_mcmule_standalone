#include "mcmule.h"
#include <string>
#include <iostream>
#include <cmath>
#include <random>
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dist(0, 100);
extern "C"
{
    // McMule will calculate this many histograms, cannot be zero but
    // very large without performance penalty.
    int mcmule_number_hist = 1;
    // The number of bins per histogram. All histograms need to have
    // the same number of bins but not all need to be useful. Can be
    // very large.
    int mcmule_number_bins = 100;
    // This defines the upper and lower bounds of the histograms
    // McMule will compute. Must be set for all histograms requested.
    double mcmule_lower_bounds[1] = {0.};
    double mcmule_upper_bounds[1] = {100.};
    // McMule will perform this many extra integrations beyond just
    // phase space etc. Useful for hit-and-miss sampling,
    // non-monochromatic beams, random acceptance etc.
    int mcmule_user_integration_dimension = 0;

    // the which_piece will either start with ee2ee (Moller) or mp2mp
    // (e-p, yes really)
    extern char __global_def_MOD_which_piece[25];

    void kelly_proton_ff(double *q2, double *Ge, double *Gm)
    {
        // This form factor cannot be used with the TPE calculation!
        double Mproton = 938.272088;
        double a11 = 2.90966;
        double a12 = -1.11542229;
        double a13 = 3.866171e-2;
        double b11 = 14.5187212;
        double b12 = 40.88333;
        double b13 = 99.999998;
        double b14 = 4.579e-5;
        double b15 = 10.3580447;

        double tau = *q2 / 4 / Mproton / Mproton;

        double tau2 = tau * tau;
        double tau3 = tau2 * tau;
        double tau4 = tau3 * tau;
        double tau5 = tau4 * tau;

        *Ge = (1 + a11 * tau + a12 * tau2 + a13 * tau3) / (1 + b11 * tau + b12 * tau2 + b13 * tau3 + b14 * tau4 + b15 * tau5);
        *Gm = 2.79284734 * *Ge;
    }

    // These are random values from fitting the Kelly form factor
    double rational11P = 2.95858e-6;
    double rational11R = 0.00445253;

    void rational11_proton_ff(double *q2, double *Ge, double *Gm)
    {
        // This form factor cannot be used with the TPE calculation!
        *Ge = (1 + (rational11P - rational11R * rational11R / 6) * *q2) / (1 + rational11P * *q2);
        *Gm = 2.79284734 * *Ge;
    }

    void mcmule_user_initialisation(void)
    {
        // This gets called on McMule start up and needs to perform
        // all initialisation, both of McMule's flavour system as well
        // as Geant4
        double Mel = 0.510998950;    // MeV
        double Mproton = 938.272088; // MeV
        double Ebeam = 1.1e3;        // MeV

        mcmule_protonff_kappa = 2.79284734;
        mcmule_protonff_lambda = 0.71e6; // Lambda^2 in MeV
                                         // std::cout << "OPENLOOPS DIR: " << __OPENLOOPS_INSTALL_DIR << std::endl;
        double s;
        std::cout << "We are running Moller scattering with Ebeam = " << Ebeam << std::endl;
        s = 2 * Mel * Mel + 2 * Mel * Ebeam;
        std::cout << "==>> Initialization of McMule...";
        mcmule_initflavour("e-p", &s);
        std::cout << "Done." << std::endl;
    }

    void boost_rf(double *rec, double *mo1)
    {
        // Boosts mo1 into the frame of rec. in-place
        double mass = rec[3] * rec[3] - rec[0] * rec[0] - rec[1] * rec[1] - rec[2] * rec[2];
        mass = std::sqrt(mass);
        double energy = rec[3];
        double dot_dot = rec[0] * mo1[0] + rec[1] * mo1[1] + rec[2] * mo1[2];

        mo1[0] += rec[0] * (dot_dot / (energy + mass) - mo1[3]) / mass;
        mo1[1] += rec[1] * (dot_dot / (energy + mass) - mo1[3]) / mass;
        mo1[2] += rec[2] * (dot_dot / (energy + mass) - mo1[3]) / mass;
        mo1[3] = (energy * mo1[3] - dot_dot) / mass;
    }

    void mcmule_measurement_function(double **res, double *p1, double *p2, double *p3, double *p4, double *p5, double *p6, double *p7)
    {
        // std::cout << "User measurement function starts" << std::endl;
        // This gets called for each event. p1, ..., p7 are pointers
        // to momenta in (px, py, pz, E). p1 and p2 are the incoming
        // particle which we will ignore. The return values are set
        // using res[0][i] to set the i-th histogram. We use the
        // mcmule_pass_cut array to accept or reject events.

        // The particles are generated in the CMS frame so we need to
        // perform a boost into the rest-frame of p2

        boost_rf(p2, p3);
        boost_rf(p2, p4);
        boost_rf(p2, p5);
        boost_rf(p2, p6);
        boost_rf(p2, p7);
        // std::cout << "Boosting to lab frame done" << std::endl;

        MCMULE_SET_NAME(0, "E");
	res[0][0] = dist(gen);
        if (res[0][0] > 0 && res[0][0] < 50)
        {
            for (int i = 0; i < mcmule_number_hist; ++i)
            {
                mcmule_pass_cut[i] = true;
            }
        }
        else
        {
            for (int i = 0; i < mcmule_number_hist; ++i)
            {
                mcmule_pass_cut[i] = false;
            }
        }
        

        
        // std::cout << "User measurement function finished" << std::endl;
        // primaryGen->ClearVector();
    }

    void mcmule_user_integration(double *x, int *ndim)
    {
        // x is a list of ndim random numbers between 0 and 1 we can
        // do with what we like.

        // mcmule_userweight can be used for partial acceptance. The
        // integral
        //    \int d\sigma -> \int d\simga * mcmule_userweight

        mcmule_userweight = 1.;
    }
};
