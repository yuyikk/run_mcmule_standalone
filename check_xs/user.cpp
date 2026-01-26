#include "mcmule.h"
#include <string>
#include <iostream>
#include <random>
#include <chrono>
int c = 0;
uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 gen(seed);                                
std::uniform_real_distribution<double> dist(0.0, 1.0); 
extern "C"
{
    // McMule will calculate this many histograms, cannot be zero but
    // very large without performance penalty.
    int mcmule_number_hist = 3;
    size_t mcmule_namelen = 10;
    // The number of bins per histogram. All histograms need to have
    // the same number of bins but not all need to be useful. Can be
    // very large.
    int mcmule_number_bins = 745;
    // This defines the upper and lower bounds of the histograms
    // McMule will compute. Must be set for all histograms requested.
    double mcmule_lower_bounds[3] = {0.55, 0.55, 0.55};
    double mcmule_upper_bounds[3] = {8.00, 8.00, 8.00};
    // McMule will perform this many extra integrations beyond just
    // phase space etc. Useful for hit-and-miss sampling,
    // non-monochromatic beams, random acceptance etc.
    int mcmule_user_integration_dimension = 0;

    // the which_piece will either start with ee2ee (Moller) or mp2mp
    // (e-p, yes really)
    extern char __global_def_MOD_which_piece[25];
    const double pi = 3.14159265358979323846;

    const double rad2deg = 180. / pi;
    int fortran_read_int();
    double fortran_read_real();
    // some helper functions
    double dot3D(const double *p1, const double *p2)
    {
        return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
    }
    double dot4D(const double *p1, const double *p2)
    {
        return p1[3] * p2[3] - dot3D(p1, p2);
    }
    double get_theta(const double *p1, const double *p2)
    {
        double dot = dot3D(p1, p2);
        double mag1 = std::sqrt(dot3D(p1, p1));
        double mag2 = std::sqrt(dot3D(p2, p2));
        if (mag1 == 0 || mag2 == 0)
        {
            std::cout << "Warning: zero-length momentum vector!" << std::endl;
            return 0.;
        }
        double cos_th = dot / mag1 / mag2;
        if (cos_th > 1)
        {
            cos_th = 1.;
        }
        else if (cos_th < -1)
        {
            cos_th = -1.;
        }
        return std::acos(cos_th);
    }
    double get_Q2(const double *p1, const double *p2);

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
    double factor = 25.68189504e-6; // conversion factor for fm^2 to MeV^-2

    void rational11_proton_ff(double *q2, double *Ge, double *Gm)
    {
        // This form factor cannot be used with the TPE calculation!
        *Ge = (1 + (rational11P - rational11R * rational11R / 6) * *q2) / (1 + rational11P * *q2);
        *Gm = 2.79284734 * *Ge;
    }

    void mcmule_user_initialisation(void)
    {

        double Mel = 0.510998950;    // MeV
        double Mproton = 938.272088; // MeV

        mcmule_protonff_kappa = 2.79284734;
        mcmule_protonff_lambda = 0.71e6; // Lambda^2 in MeV
        std::cout << "Enter beam energy(MeV): " << std::endl;
        double Ebeam = fortran_read_real();
        double s;
        if (__global_def_MOD_which_piece[1] == 'e')
        {
            std::cout << "We are running Moller scattering with Ebeam (MeV) = " << Ebeam << std::endl;
            s = 2 * Mel * Mel + 2 * Mel * Ebeam;
            mcmule_number_hist = 2;
        }
        else if (__global_def_MOD_which_piece[1] == 'p')
        {
            std::cout << "We are running e-p scattering with Ebeam (MeV) = " << Ebeam << std::endl;
            s = Mel * Mel + Mproton * Mproton + 2 * Mproton * Ebeam;
            std::cout << "Please choose proton form factor models:\n";
            std::cout << "  0: dipole form factor.\n";
            std::cout << "  1: monopole form factor.\n";
            std::cout << "  2: Kelly form factor.\n";
            std::cout << "  3: Rational(1, 1) form factor.\n";
            int user_input = fortran_read_int();
            // std::cin >> user_input;
            if (user_input == 0)
            {
                std::cout << "Using McMule's default dipole form factors" << std::endl;
                std::cout << "with lambda = " << mcmule_protonff_lambda << " and kappa = " << mcmule_protonff_kappa << std::endl;
                mcmule_protonff = &__nucl_protonff_MOD_proton_dipole;
            }
            else if (user_input == 1)
            {
                std::cout << "Using McMule's default monopole form factors" << std::endl;
                std::cout << "with lambda = " << mcmule_protonff_lambda << " and kappa = " << mcmule_protonff_kappa << std::endl;
                mcmule_protonff = &__nucl_protonff_MOD_proton_monopole;
            }
            else if (user_input == 2)
            {
                std::cout << "Using Kelly form factor" << std::endl;
                mcmule_protonff = &kelly_proton_ff;
            }
            else if (user_input == 3)
            {
                std::cout << "Please enter parameter R(adius): " << std::endl;
                rational11R = fortran_read_real() * factor;
                std::cout << "Please enter parameter P: " << std::endl;
                rational11P = fortran_read_real() * factor;
                std::cout << "Using Rational(1,1) form factor" << std::endl;
                std::cout << "with p = " << rational11P / factor << " and R = " << rational11R / factor << std::endl;
                mcmule_protonff = &rational11_proton_ff;
            }
        }
        if (!mcmule_protonff)
        {
            std::cerr << "Invalid form factor" << std::endl;
        }

        mcmule_initflavour("e-p", &s);
    }

    void boost_rf(double *rec, double *mo1)
    {
        // Boosts mo1 into the frame of rec. in-place
        double mass = rec[3] * rec[3] - rec[0] * rec[0] - rec[1] * rec[1] - rec[2] * rec[2];
        mass = sqrt(mass);
        double energy = rec[3];
        double dot_dot = rec[0] * mo1[0] + rec[1] * mo1[1] + rec[2] * mo1[2];

        mo1[0] += rec[0] * (dot_dot / (energy + mass) - mo1[3]) / mass;
        mo1[1] += rec[1] * (dot_dot / (energy + mass) - mo1[3]) / mass;
        mo1[2] += rec[2] * (dot_dot / (energy + mass) - mo1[3]) / mass;
        mo1[3] = (energy * mo1[3] - dot_dot) / mass;
    }

    void mcmule_measurement_function(double **res, double *p1, double *p2, double *p3, double *p4, double *p5, double *p6, double *p7)
    {
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
        double z_unit[4] = {0, 0, 1, 1};
        double theta_mc = get_theta(z_unit, p3) * rad2deg;
        double rand = dist(gen);

        if (theta_mc < 0.5 || theta_mc > 8.5)
        {
            for (int i = 0; i < mcmule_number_hist; ++i)
            {
                mcmule_pass_cut[i] = false;
            }
            return;
        }
        MCMULE_SET_NAME(0, "th_ep_01");
        MCMULE_SET_NAME(1, "th_ep_02");
        MCMULE_SET_NAME(2, "th_ep_03");
        res[0][1] = theta_mc; // this line is duplicate
        res[0][2] = theta_mc; // this line is duplicate
        res[0][3] = theta_mc; // this line is duplicate
        mcmule_pass_cut[0] = (theta_mc > 0.5 && theta_mc < 8.5);
        mcmule_pass_cut[1] = (mcmule_pass_cut[1] && rand < 0.7);
        mcmule_pass_cut[2] = (mcmule_pass_cut[1] && rand < 0.3);


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
