
#include "mcmule.h"
#include <iostream>
#include <string>
#include <cmath>
extern "C"
{
    // McMule will calculate this many histograms, cannot be zero but
    // very large without performance penalty.
    const int mcmule_number_hist = 12;
    // The number of bins per histogram. All histograms need to have
    // the same number of bins but not all need to be useful. Can be
    // very large.
    const int mcmule_number_bins = 250;
    // This defines the upper and lower bounds of the histograms
    // McMule will compute. Must be set for all histograms requested.
    double mcmule_lower_bounds[mcmule_number_hist] = {0.5, 0.5, 0.5, 0.5,
                                                      0.0, 0.0, 0.0, 0.0,
                                                      0.5, 0.5, 0.0, 0.0};
    double mcmule_upper_bounds[mcmule_number_hist] = {6.75, 6.75, 6.75, 6.75,
                                                      2500., 2500., 2500., 2500.,
                                                      6.75, 6.75, 100.e3, 100.e3};
    // McMule will perform this many extra integrations beyond just
    // phase space etc. Useful for hit-and-miss sampling,
    // non-monochromatic beams, random acceptance etc.
    int mcmule_user_integration_dimension = 0;

    // the which_piece will either start with ee2ee (Moller) or mp2mp
    // (e-p, yes really)
    extern char __global_def_MOD_which_piece[25];
    int process = 0;
    constexpr double PI = 3.14159265358979323846;
    constexpr double deg = PI / 180.;

    double Q2(const double *p1, const double *p2)
    {
        double qx = p1[0] - p2[0];
        double qy = p1[1] - p2[1];
        double qz = p1[2] - p2[2];
        double qE = p1[3] - p2[3];

        // q^2 = E^2 - |p|^2
        double q2 = qE * qE - (qx * qx + qy * qy + qz * qz);

        // Q^2 = -q^2
        return -q2;
    }
    double angle_between(const double *a, const double *b)
    {
        // dot product of spatial components
        double dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

        // magnitudes of spatial vectors
        double magA = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        double magB = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);

        if (magA == 0.0 || magB == 0.0)
            return 0.0; // avoid division by zero

        // cosine of the angle
        double cosTheta = dot / (magA * magB);

        // numerical safeguard (due to floating-point precision)
        if (cosTheta > 1.0)
            cosTheta = 1.0;
        if (cosTheta < -1.0)
            cosTheta = -1.0;

        // return angle in radians
        return std::acos(cosTheta);
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
        std::string which_piece(__global_def_MOD_which_piece);
        std::string which_piece_sub = which_piece.substr(0, 5);
	std::cout << "which_piece_sub: " << which_piece_sub << std::endl;
        double s;
        if (process == 0)
        {
            std::cout << "==> Running e-p elastic scattering with Ebeam " << Ebeam << " MeV." << std::endl;
            std::cout << "      additional cut on photons with E_y > 20 MeV if th(el_f,y)>6 mrad" << std::endl;
            if (which_piece_sub == "mp2mp")
            {
                s = Mel * Mel + Mproton * Mproton + 2 * Mproton * Ebeam;
            }
        }
        else if (process == 1)
        {
            std::cout << "==> Running Moller scattering with Ebeam " << Ebeam << " MeV." << std::endl;
            if (which_piece_sub == "ee2ee")
            {
                s = 2 * Mel * Mel + 2 * Mel * Ebeam;
            }
        }
        else
        {
            std::cout << "Invalid process code." << std::endl;
            return;
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

        double p1_rest[4] = {p1[0], p1[1], p1[2], p1[3]};
        double p2_rest[4] = {p2[0], p2[1], p2[2], p2[3]};
        double p3_rest[4] = {p3[0], p3[1], p3[2], p3[3]};
        double p4_rest[4] = {p4[0], p4[1], p4[2], p4[3]};
        double p5_rest[4] = {p5[0], p5[1], p5[2], p5[3]};
        double p6_rest[4] = {p6[0], p6[1], p6[2], p6[3]};
        double p7_rest[4] = {p7[0], p7[1], p7[2], p7[3]};
        double qsql, qsqp;
        double thetal, cthetal, thetalJ, cthetalJ;
        double thetal1, cthetal1, thetal2, cthetal2;
        double thetay35, cthetay35, thetay36, cthetay36;
        double cone = 6e-3, cone_lim = 20., eytot;
        double th_n, th_w, th_s, th_h, e_s, e_h, e_n, e_w;


        boost_rf(p2, p1_rest);
        boost_rf(p2, p2_rest);
        boost_rf(p2, p3_rest);
        boost_rf(p2, p4_rest);
        boost_rf(p2, p5_rest);
        boost_rf(p2, p6_rest);
        boost_rf(p2, p7_rest);

        // Now calculate variables
        thetay35 = angle_between(p3_rest, p5_rest);
        cthetay35 = std::cos(thetay35);

        thetay36 = angle_between(p3_rest, p6_rest);
        cthetay36 = std::cos(thetay36);

        double p3_rest_j[4] = {p3_rest[0], p3_rest[1], p3_rest[2], p3_rest[3]};

        if (thetay35 < cone)
        {
            for (int idx = 0; idx < 4; ++idx)
            {
                p3_rest_j[idx] = p3_rest[idx] + p5_rest[idx];
            }
        }

        if (thetay36 < cone)
        {
            for (int idx = 0; idx < 4; ++idx)
            {
                p3_rest_j[idx] = p3_rest_j[idx] + p6_rest[idx];
            }
        }

        thetal = angle_between(p1_rest, p3_rest);
        cthetal = std::cos(thetal);
        thetal1 = thetal;
        cthetal1 = cthetal;
        thetal2 = angle_between(p1_rest, p4_rest);
        cthetal2 = std::cos(thetal2);

        thetalJ = angle_between(p1_rest, p3_rest_j);
        cthetalJ = std::cos(thetalJ);
        bool pass_cut = true;
        if (process == 0)
        {
            if (thetal < 0.7 * deg || thetal > 6 * deg)
            {
                pass_cut = false;
            }
            eytot = 0;
            if (thetay35 > cone)
            {
                eytot += p5_rest[3];
            }
            if (thetay36 > cone)
            {
                eytot += p6_rest[3];
            }
            if (eytot > cone_lim)
            {
                pass_cut = false;
            }
            qsql = Q2(p1, p3);
            qsqp = Q2(p2, p4);
            MCMULE_SET_NAME(8, "th_l");
            res[0][8] = thetal / deg;
            MCMULE_SET_NAME(9, "th_lJ");
            res[0][9] = thetalJ / deg;
            MCMULE_SET_NAME(10, "qsql");
            res[0][10] = qsql;
            MCMULE_SET_NAME(11, "qsqp");
            res[0][11] = qsqp;
        }
        else if (process == 1)
        {
            if (thetal1 < 0.7 * deg || thetal1 > 6 * deg || thetal2 < 0.7 * deg || thetal2 > 6 * deg)
            {
                pass_cut = false;
            }
            double th3, th4;
            // TVector3 unit_z(0, 0, 1);
            double unit_z[4] = {0, 0, 1, 0};
            th3 = angle_between(p3_rest, unit_z);
            th4 = angle_between(p4_rest, unit_z);
            double e3 = p3_rest[3];
            double e4 = p4_rest[3];
            if (th3 <= th4)
            {
                th_n = th3;
                th_w = th4;
                e_n = e3;
                e_w = e4;
            }
            else
            {
                th_n = th4;
                th_w = th3;
                e_n = e4;
                e_w = e3;
            }
            if (e3 <= e4)
            {
                th_s = th3;
                e_s = e3;
                th_h = th4;
                e_h = e4;
            }
            else
            {
                th_s = th4;
                e_s = e4;
                th_h = th3;
                e_h = e3;
            }

            MCMULE_SET_NAME(0, "th_n");
            res[0][0] = th_n;
            MCMULE_SET_NAME(1, "th_w");
            res[0][1] = th_w;
            MCMULE_SET_NAME(2, "th_s");
            res[0][2] = th_s;
            MCMULE_SET_NAME(3, "th_h");
            res[0][3] = th_h;
            MCMULE_SET_NAME(4, "e_n");
            res[0][4] = e_n;
            MCMULE_SET_NAME(5, "e_w");
            res[0][5] = e_w;
            MCMULE_SET_NAME(6, "e_s");
            res[0][6] = e_s;
            MCMULE_SET_NAME(7, "e_h");
            res[0][7] = e_h;
        }
	    for (int i = 0; i < mcmule_number_hist; ++i)
            {
             	mcmule_pass_cut[i] = pass_cut;
            }
    }

    void mcmule_user_integration(double *x, int *ndim)
    {
        // x is a list of ndim random numbers between 0 and 1 we can
        // do with what we like.

        // mcmule_userweight can be used for partial acceptance. The
        // integral
        //    \int d\sigma . \int d\simga * mcmule_userweight

        mcmule_userweight = 1.;
    }
};
