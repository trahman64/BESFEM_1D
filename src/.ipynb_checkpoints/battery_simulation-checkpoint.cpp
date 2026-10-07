#include "mfem.hpp"

#include "../includes/spherical_diffusion.hpp"
#include "../includes/linear_diffusion.hpp"
#include "../includes/stat_potential.hpp"
#include "../includes/diffCoeff_utils.hpp"
#include "../includes/atnv_calculator.hpp"
#include "../includes/reaction.hpp"

#include <iostream>
#include <vector>


// ============================================================
// Constants
// ============================================================

const double F = 96485.332;
const double R = 8.314;
const double T = 300.0;

const double alpha = 0.5;

const double i0 = 0.5e-3;

const double a = 2.409e3;
const double t_minus = 0.7619;

const double rho = 0.0312;


// Volume fraction of liquid
const double eps_l_eld = 0.301;
const double eps_l_sep = 1.0;

// Tortuosity of liquid
const double tau_l_eld = 1.521;
const double tau_l_sep = 1.0;

const double C0 = 0.001;
const double De = 0.25e-5;


// Electrolyte time step
const double dt = 1.0e-4;


// Solid potential
double kappa_s = 0.075;

// Volume fraction of solid
const double eps_s_eld = 0.699;
const double eps_s_sep = 1.0e0;

// Tortuosity of solid
const double tau_s_eld = 1.324;
const double tau_s_sep = 1.0e-3;


// Solid potential boundary value
double BvP = 0.0;


// ============================================================
// Main
// ============================================================

int main()
{
    // ========================================================
    // Mesh and finite element space
    // ========================================================

    mfem::Mesh mesh("../inputs/Mesh_80_1D02.mesh");

    int dim = mesh.Dimension();
    int order = 1;

    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fespace(&mesh, &fec);


    // ========================================================
    // Electrolyte diffusion
    // ========================================================

    mfem::GridFunction source_ely(&fespace);
    source_ely = 0.0;

    std::cout << "Creating linear diffusion" << std::endl;

    LinearDiffusion salt_electrolyte(
        &mesh,
        &fespace,
        eps_l_sep,
        eps_l_eld,
        De,
        tau_l_sep,
        tau_l_eld,
        C0,
        dt
    );


    // Current electrolyte concentration
    mfem::GridFunction Ce(&fespace);
    Ce = salt_electrolyte.GetConcentration();


    // ========================================================
    // Electrolyte diffusion coefficient
    // ========================================================

    mfem::GridFunction De_gf(&fespace);

    De_gf = ComputeDamb(Ce);


    // ========================================================
    // Spherical particles
    // ========================================================

    std::cout << "Creating spherical particles" << std::endl;

    int n = 60;

    std::vector<SphericalDiffusion> particles;

    particles.reserve(n);

    for (int i = 0; i < n; i++)
    {
        particles.emplace_back(
            4.0e-4,       // particle radius
            40,           // radial elements
            1.5e-10,      // diffusion coefficient
            i,            // particle ID
            1.0e-2,       // particle dt
            1,            // FE order
            0.3            // initial concentration
        );
    }


    // Surface concentration field used by Butler-Volmer
    mfem::GridFunction Cp_surf(&fespace);

    Cp_surf = 0.0;


    // Initialize particle surface concentrations
    for (int j = 0; j < n; j++)
    {
        int index = 20 + j;

        Cp_surf(index) =
            particles[j].GetMeanConcentration();
    }


    // ========================================================
    // Solid potential
    // ========================================================

    std::cout << "Creating solid potential" << std::endl;

    mfem::Array<int> ess_bdr_s(mesh.bdr_attributes.Max());

    ess_bdr_s = 0;

    // Right boundary
    ess_bdr_s[1] = 1;


    mfem::GridFunction Kappa(&fespace);

    Kappa = kappa_s;


    StatPotential solid_potential(
        &fespace,
        ess_bdr_s
    );


    solid_potential.SetWeightVector(
        eps_s_sep,
        eps_s_eld,
        tau_s_sep,
        tau_s_eld
    );


    solid_potential.SetCoefficient(Kappa);

    solid_potential.BuildOperator();


    // Additional solid potential term
    mfem::Vector AtnV_0(
        fespace.GetTrueVSize()
    );

    AtnV_0 = 0.0;


    // Solid potential source
    mfem::GridFunction source_phs(&fespace);

    source_phs = 0.0;


    // ========================================================
    // Liquid potential
    // ========================================================

    std::cout << "Creating liquid potential" << std::endl;


    mfem::Array<int> ess_bdr_l(
        mesh.bdr_attributes.Max()
    );

    ess_bdr_l = 0;

    // Left boundary
    ess_bdr_l[0] = 1;


    double Cst1 =
        F / R / T;


    double tc1 =
        (2.0 * t_minus - 1.0)
        /
        (2.0 * t_minus * (1.0 - t_minus));


    double tc2 =
        1.0
        /
        (2.0 * t_minus * (1.0 - t_minus))
        *
        Cst1;


    double scaleConst =
        tc2 / tc1 * Cst1;


    mfem::GridFunction Dmp(&fespace);

    mfem::GridFunction Kpl(&fespace);


    Dmp = De_gf;

    Dmp *= tc1;


    Kpl = De_gf;

    Kpl *= scaleConst;

    Kpl *= Ce;


    StatPotential liquid_potential(
        &fespace,
        ess_bdr_l
    );


    liquid_potential.SetWeightVector(
        eps_l_sep,
        eps_l_eld,
        tau_l_sep,
        tau_l_eld
    );


    liquid_potential.SetCoefficient(Kpl);

    liquid_potential.BuildOperator();


    // ========================================================
    // References to potentials
    //
    // IMPORTANT:
    // These are references to the GridFunctions stored inside
    // StatPotential. When Solve() changes phi, these references
    // automatically see the new values.
    // ========================================================

    mfem::GridFunction &phi_p =
        solid_potential.GetPotential();


    mfem::GridFunction &phi_e =
        liquid_potential.GetPotential();


    // ========================================================
    // AtnV term for liquid potential
    // ========================================================

    AtnVCalculator AtnVCalt(
        &fespace,
        Dmp
    );


    mfem::Vector &AtnV =
        AtnVCalt.Compute(Ce);


    // Liquid potential source
    mfem::GridFunction source_phl(&fespace);

    source_phl = 0.0;


    double BvE = 0.0;


    // ========================================================
    // Butler-Volmer
    // ========================================================

    std::cout << "Creating Butler-Volmer" << std::endl;


    ButlerVolmer bv(
        &fespace,
        i0,
        a,
        alpha,
        F,
        R,
        T
    );


    // Reference to reaction stored inside ButlerVolmer.
    //
    // When bv.Compute() updates its rxn GridFunction,
    // this reference automatically sees the new values.

    mfem::GridFunction &rxn =
        bv.GetReaction();


    // ========================================================
    // INITIAL POTENTIAL SOLVE
    //
    // Butler-Volmer requires phi_p and phi_e.
    //
    // However, phi_p and phi_e are initially zero because
    // StatPotential initializes phi = 0.
    //
    // Therefore, use the old prescribed reaction ONCE to
    // calculate initial potentials.
    // ========================================================

    std::cout << "Calculating initial potentials" << std::endl;


    mfem::GridFunction initial_rxn(&fespace);

    initial_rxn = 0.0;


    // 60 electrode locations:
    // 20, 21, ..., 79

    for (int i = 20; i < 80; i++)
    {
        initial_rxn(i) =
            a * 0.02e-6;
    }


    // --------------------------------------------------------
    // Initial solid potential
    // --------------------------------------------------------

    source_phs = initial_rxn;

    source_phs *= F;


    solid_potential.Solve(
        source_phs,
        AtnV_0,
        BvP
    );


    // --------------------------------------------------------
    // Initial liquid potential
    // --------------------------------------------------------

    source_phl = initial_rxn;


    liquid_potential.Solve(
        source_phl,
        AtnV,
        BvE
    );


    // phi_p and phi_e now contain the solutions because they
    // are references to the GridFunctions inside StatPotential.


    // ========================================================
    // Time loop
    // ========================================================

    std::cout << "Starting loop" << std::endl;


    int num_steps = 10000;


    for (int step = 0; step < num_steps; step++)
    {
        // ====================================================
        // 1. CURRENT ELECTROLYTE CONCENTRATION
        // ====================================================

        Ce = salt_electrolyte.GetConcentration();


        // ====================================================
        // 2. CURRENT PARTICLE SURFACE CONCENTRATIONS
        // ====================================================

        Cp_surf = 0.0;


        for (int j = 0; j < n; j++)
        {
            int index = 20 + j;


            Cp_surf(index) =
                particles[j].GetMeanConcentration();
        }


        // ====================================================
        // 3. BUTLER-VOLMER
        //
        // Inputs:
        //     Ce
        //     Cp_surf
        //     phi_e
        //     phi_p
        //
        // Output:
        //     rxn
        // ====================================================

        bv.Compute(
            Ce,
            Cp_surf,
            phi_e,
            phi_p
        );


        // rxn now contains the NEW reaction calculated from
        // Butler-Volmer.


        // ====================================================
        // 4. BUILD SOURCES FROM NEW REACTION
        // ====================================================


        // ----------------------------------------------------
        // Electrolyte concentration source
        // ----------------------------------------------------

        source_ely = rxn;

        source_ely *= -t_minus;


        // ----------------------------------------------------
        // Solid potential source
        // ----------------------------------------------------

        source_phs = rxn;

        source_phs *= F;


        // ----------------------------------------------------
        // Liquid potential source
        // ----------------------------------------------------

        source_phl = rxn;


        // ====================================================
        // 5. UPDATE ELECTROLYTE CONCENTRATION
        // ====================================================

        salt_electrolyte.Stepping(
            source_ely
        );


        // ====================================================
        // 6. UPDATE SOLID POTENTIAL
        // ====================================================

        solid_potential.Solve(
            source_phs,
            AtnV_0,
            BvP
        );


        // phi_p is automatically updated because it is a
        // reference to solid_potential's internal phi.


        // ====================================================
        // 7. UPDATE LIQUID POTENTIAL
        // ====================================================

        liquid_potential.Solve(
            source_phl,
            AtnV,
            BvE
        );


        // phi_e is automatically updated because it is a
        // reference to liquid_potential's internal phi.


        // ====================================================
        // 8. UPDATE SPHERICAL PARTICLES
        // ====================================================

        for (int j = 0; j < n; j++)
        {
            int index = 20 + j;


            double surface_flux =
                rxn(index) / a;


            particles[j].Stepping(
                surface_flux
            );
        }


        // ====================================================
        // Optional debugging
        // ====================================================

        if (step % 100 == 0)
        {
            std::cout
                << "Step "
                << step
                << std::endl;

            std::cout
                << "rxn(20) = "
                << rxn(20)
                << std::endl;

            std::cout
                << "Cp_surf(20) = "
                << Cp_surf(20)
                << std::endl;

            std::cout
                << "phi_p(20) = "
                << phi_p(20)
                << std::endl;

            std::cout
                << "phi_e(20) = "
                << phi_e(20)
                << std::endl;
        }
    }


    // ========================================================
    // Save / print results
    // ========================================================

    Ce = salt_electrolyte.GetConcentration();


    Ce.Print();


    salt_electrolyte.SaveConc();


    double MnConc =
        salt_electrolyte.GetMeanConcentration();


    std::cout
        << "MC = "
        << MnConc
        << std::endl;


    liquid_potential.SavePote();


    return 0;
}