/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dfChemistryModel.H"
#include "UniformField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::dfChemistryModel<ThermoType>::dfChemistryModel
(
    ThermoType& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            "CanteraTorchProperties",
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    thermo_(thermo),
    mixture_(dynamic_cast<CanteraMixture&>(thermo)),
    CanteraGas_(mixture_.CanteraGas()),
    mesh_(thermo.p().mesh()),
    chemistry_(lookup("chemistry")),
    relTol_(this->subDict("odeCoeffs").lookupOrDefault("relTol",1e-9)),
    absTol_(this->subDict("odeCoeffs").lookupOrDefault("absTol",1e-15)),
    Y_(mixture_.Y()),
    rhoD_(mixture_.nSpecies()),
    hai_(mixture_.nSpecies()),
    yTemp_(mixture_.nSpecies()),
    dTemp_(mixture_.nSpecies()),
    hrtTemp_(mixture_.nSpecies()),
    cTemp_(mixture_.nSpecies()),
    RR_(mixture_.nSpecies()),
    alpha_(const_cast<volScalarField&>(thermo.alpha())),
    T_(thermo.T()),
    p_(thermo.p()),
    rho_(mesh_.objectRegistry::lookupObject<volScalarField>("rho")),
    mu_(const_cast<volScalarField&>(dynamic_cast<psiThermo&>(thermo).mu()())),
    psi_(const_cast<volScalarField&>(dynamic_cast<psiThermo&>(thermo).psi())),
    Qdot_
    (
        IOobject
        (
            "Qdot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    ),
    torchSwitch_(lookupOrDefault("torch", false))
{
    if(torchSwitch_)
    {
        torchModelName_ = this->lookupOrDefault("torchModel", word(""));
        Xmu_ = scalarList(this->subDict("torchParameters").lookup("Xmu"));
        Xstd_ = scalarList(this->subDict("torchParameters").lookup("Xstd"));
        Ymu_ = scalarList(this->subDict("torchParameters").lookup("Ymu"));
        Ystd_ = scalarList(this->subDict("torchParameters").lookup("Ystd"));
        Tact_ = this->subDict("torchParameters").lookupOrDefault("Tact", 700);
        Qdotact_ = this->subDict("torchParameters").lookupOrDefault("Qdotact", 1e9);
    }

    for(const auto& name : CanteraGas_->speciesNames())
    {
        species_.append(name);
    }
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    forAll(rhoD_, i)
    {
        rhoD_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoD_" + Y_[i].name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimDensity*dimViscosity, 0) // kg/m/s
            )
        );
    }

    forAll(hai_, i)
    {
        hai_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hai_" + Y_[i].name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimEnergy/dimMass, 0)
            )
        );
    }

    Info<<"--- I am here in Cantera-construct ---"<<endl;
    Info<<"relTol_ === "<<relTol_<<endl;
    Info<<"absTol_ === "<<absTol_<<endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::dfChemistryModel<ThermoType>::
~dfChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::dfChemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    scalar result = 0;
    if(torchSwitch_)
    {
        result = torchSolve(deltaT);
    }
    else
    {
        result = canteraSolve(deltaT);
    }
    return result;
}

template<class ThermoType>
Foam::scalar Foam::dfChemistryModel<ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ThermoType>
Foam::scalar Foam::dfChemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::dfChemistryModel<ThermoType>::canteraSolve
(
    const DeltaTType& deltaT
)
{
    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    Info<<"=== begin Cantera-solve === "<<endl;

    Cantera::Reactor react;
    //Cantera::IdealGasReactor react;  // Constant-UV, default, constant Volumn
    //Cantera::IdealGasConstPressureReactor react;  // Constant-HP, constant pressure

    scalarField c0(CanteraGas_->nSpecies());

    Qdot_ = Zero;

    forAll(T_, cellI)
    {
        scalar Ti = T_[cellI];
        scalar pi = p_[cellI];
        try
        {
            for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
            {
                yTemp_[i] = Y_[i][cellI];
            }

            CanteraGas_->setState_TPY(Ti, pi, yTemp_.begin());
            CanteraGas_->getConcentrations(c0.begin()); // value --> c0

            react.insert(mixture_.CanteraSolution());
            react.setEnergy(0); // keep T const before and after sim.advance. this will give you a little improvement
            Cantera::ReactorNet sim;
            sim.addReactor(react);
            setNumerics(sim);


            sim.advance(deltaT[cellI]);


            CanteraGas_->getConcentrations(cTemp_.begin()); // value --> cTemp_

            for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
            {
                RR_[i][cellI] = (cTemp_[i] - c0[i])*CanteraGas_->molecularWeight(i)/deltaT[cellI];
            }
            // CanteraGas_->molecularWeight(i)    kg/kmol

            forAll(Y_, i)
            {
                const scalar hc = CanteraGas_->Hf298SS(i)/CanteraGas_->molecularWeight(i); // J/kg
                Qdot_[cellI] -= hc*RR_[i][cellI];
            }
        }
        catch(Cantera::CanteraError& err)
        {
            // handle exceptions thrown by Cantera
            std::cout << err.what() << std::endl;

            FatalErrorIn("dfChemistryModel::solve")
                << " Cantera complained in cell " << cellI
                << " with a Cantera::CanteraError"  << endl
                << abort(FatalError) ;
        }
    }

    Info<<"=== end Cantera-solve === "<<endl;
    return deltaTMin;
}


template<class ThermoType>
void Foam::dfChemistryModel<ThermoType>::setNumerics(Cantera::ReactorNet &sim)
{
    sim.setTolerances(relTol_,absTol_);
}


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::dfChemistryModel<ThermoType>::torchSolve
(
    const DeltaTType& deltaT
)
{
    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    Info<<"=== begin torch&ode-solve === "<<endl;

    // set variables
    scalarList yPre_(mixture_.nSpecies());
    scalarList yBCT_(mixture_.nSpecies());
    scalarList u_(mixture_.nSpecies()+2); //plus T and p
    Cantera::Reactor react;
    double lambda = 0.1;

    // load model
    torch::DeviceType device_type;
    torch::DeviceIndex index = -1;
    if (0)
    {
        device_type = torch::kCUDA;
        index = 0;
    }
    else
    {
        device_type = torch::kCPU;
        index = 0;
    }
    torch::Device device(device_type, index);
    if (!isFile(torchModelName_))
    {
        FatalErrorInFunction
            << torchModelName_
            << " doesn't exist!"
            << exit(FatalError);
    }
    torch::jit::script::Module torchModel_ = torch::jit::load(torchModelName_, device);

    std::vector<size_t> torch_cell;
    label torch_cellname= 0;

    // obtain the number of DNN cells
    forAll(T_, cellI)
    {
        if ((Qdot_[cellI] >= Qdotact_) && (T_[cellI] >= Tact_))
        {
            torch_cell.push_back(cellI);
        }
    }

    torch::Tensor inputs(torch::zeros({torch_cell.size(),mixture_.nSpecies()+2}));
    forAll(T_, cellI)
    {
        scalar Ti = T_[cellI];
        scalar pi = p_[cellI];
        scalar rhoi = rho_[cellI];

        if ((Qdot_[cellI] >= Qdotact_) && (Ti >= Tact_))
        {
            Qdot_[cellI] = 0.0;
            // Info<<"Now is DNN "<<endl;

            // set inputs
            std::vector<double> inputs_;
            inputs_.push_back((Ti - Xmu_[0])/Xstd_[0]);
            inputs_.push_back((pi / 101325 - Xmu_[1])/Xstd_[1]);
            for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
            {
                yPre_[i] = Y_[i][cellI];
                yBCT_[i] = (pow(yPre_[i],lambda) - 1) / lambda; // function BCT
            }
            for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
            {
                inputs_.push_back((yBCT_[i] - Xmu_[i+2]) / Xstd_[i+2]);
            }
            inputs[torch_cellname] = torch::tensor(inputs_);
            torch_cellname ++;
            // Info << "cell_name = "<< cellI <<endl;
        }
        else
        {
            Qdot_[cellI] = 0.0;
            for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
            {
                yPre_[i] = Y_[i][cellI];
            }

            CanteraGas_->setState_TPY(Ti, pi, yPre_.begin());
            react.insert(mixture_.CanteraSolution());
            react.setEnergy(0);

            Cantera::ReactorNet sim;
            sim.addReactor(react);
            setNumerics(sim);
            sim.advance(deltaT);

            CanteraGas_->getMassFractions(yTemp_.begin());

            for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
            {
                RR_[i][cellI] = (yTemp_[i] - yPre_[i])*rhoi/deltaT;
                const scalar hc = CanteraGas_->Hf298SS(i)/CanteraGas_->molecularWeight(i); // J/kg
                Qdot_[cellI] -= hc*RR_[i][cellI];
            }
        }
    }
    // DNN
    std::vector<torch::jit::IValue> INPUTS{inputs};
    auto outputs_raw = torchModel_.forward(INPUTS);
    auto outputs = outputs_raw.toTensor();

    for(size_t cellI = 0; cellI<torch_cell.size(); cellI ++)
    {
        // update y
        scalar Yt = 0;
        for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
        {
            yPre_[i] = Y_[i][torch_cell[cellI]];
            yTemp_[i] = Y_[i][torch_cell[cellI]];
            yBCT_[i] = (pow(yPre_[i],lambda) - 1) / lambda; // function BCT
        }
        for (size_t i=0; i<(CanteraGas_->nSpecies()); i++)//
        {
            u_[i+2] = outputs[cellI][i+2].item().to<double>()*Ystd_[i+2]+Ymu_[i+2];
            yTemp_[i] = pow((yBCT_[i] + u_[i+2]*deltaT)*lambda+1,1/lambda);
            Yt += yTemp_[i];
        }
        for (size_t i=0; i<CanteraGas_->nSpecies(); i++)
        {
            yTemp_[i] = yTemp_[i] / Yt;
            RR_[i][torch_cell[cellI]] = (yTemp_[i] - Y_[i][torch_cell[cellI]])*rho_[torch_cell[cellI]]/deltaT;
            const scalar hc = CanteraGas_->Hf298SS(i)/CanteraGas_->molecularWeight(i); // J/kg
            Qdot_[torch_cell[cellI]] -= hc*RR_[i][torch_cell[cellI]];
        }
    }

    Info<<"=== end torch&ode-solve === "<<endl;
    return deltaTMin;
}

template<class ThermoType>
void Foam::dfChemistryModel<ThermoType>::correctThermo()
{
    forAll(T_, celli)
    {
        forAll(Y_, i)
        {
            yTemp_[i] = Y_[i][celli];
        }
        CanteraGas_->setState_PY(p_[celli], yTemp_.begin());
        CanteraGas_->setState_HP(thermo_.he()[celli], p_[celli]); // setState_HP needs (J/kg)

        T_[celli] = CanteraGas_->temperature();

        psi_[celli] = CanteraGas_->meanMolecularWeight()/CanteraGas_->RT(); // meanMolecularWeight() kg/kmol    RT() Joules/kmol

        mu_[celli] = mixture_.CanteraTransport()->viscosity(); // Pa-s

        alpha_[celli] = mixture_.CanteraTransport()->thermalConductivity()/(CanteraGas_->cp_mass()); // kg/(m*s)
        // thermalConductivity() W/m/K
        // cp_mass()   J/kg/K


        if (mixture_.transportModelName() == "UnityLewis")
        {
            forAll(rhoD_, i)
            {
                rhoD_[i][celli] = alpha_[celli];
            }
        }
        else
        {
            mixture_.CanteraTransport()->getMixDiffCoeffsMass(dTemp_.begin()); // m2/s

            CanteraGas_->getEnthalpy_RT(hrtTemp_.begin()); //hrtTemp_=m_h0_RT non-dimension
            // constant::physicoChemical::R.value()   J/(mol·k)
            const scalar RT = constant::physicoChemical::R.value()*1e3*T_[celli]; // J/kmol/K
            forAll(rhoD_, i)
            {
                rhoD_[i][celli] = rho_[celli]*dTemp_[i];

                // CanteraGas_->molecularWeight(i)    kg/kmol
                hai_[i][celli] = hrtTemp_[i]*RT/CanteraGas_->molecularWeight(i);
            }
        }
    }


    const volScalarField::Boundary& pBf = p_.boundaryField();

    const volScalarField::Boundary& rhoBf = rho_.boundaryField();

    volScalarField::Boundary& TBf = T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf = psi_.boundaryFieldRef();

    volScalarField::Boundary& hBf = thermo_.he().boundaryFieldRef();

    volScalarField::Boundary& muBf = mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf = alpha_.boundaryFieldRef();

    forAll(T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& ph = hBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                forAll(Y_, i)
                {
                    yTemp_[i] = Y_[i].boundaryField()[patchi][facei];
                }
                CanteraGas_->setState_TPY(pT[facei], pp[facei], yTemp_.begin());

                ph[facei] = CanteraGas_->enthalpy_mass();

                ppsi[facei] = CanteraGas_->meanMolecularWeight()/CanteraGas_->RT();

                pmu[facei] = mixture_.CanteraTransport()->viscosity();

                palpha[facei] = mixture_.CanteraTransport()->thermalConductivity()/(CanteraGas_->cp_mass());
                if (mixture_.transportModelName() == "UnityLewis")
                {
                    forAll(rhoD_, i)
                    {
                        rhoD_[i].boundaryFieldRef()[patchi][facei] = palpha[facei];
                    }
                }
                else
                {
                    mixture_.CanteraTransport()->getMixDiffCoeffsMass(dTemp_.begin());

                    CanteraGas_->getEnthalpy_RT(hrtTemp_.begin());
                    const scalar RT = constant::physicoChemical::R.value()*1e3*pT[facei];
                    forAll(rhoD_, i)
                    {
                        rhoD_[i].boundaryFieldRef()[patchi][facei] = prho[facei]*dTemp_[i];

                        hai_[i].boundaryFieldRef()[patchi][facei] = hrtTemp_[i]*RT/CanteraGas_->molecularWeight(i);
                    }
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                forAll(Y_, i)
                {
                    yTemp_[i] = Y_[i].boundaryField()[patchi][facei];
                }
                CanteraGas_->setState_PY(pp[facei], yTemp_.begin());
                CanteraGas_->setState_HP(ph[facei], pp[facei]);

                pT[facei] = CanteraGas_->temperature();

                ppsi[facei] = CanteraGas_->meanMolecularWeight()/CanteraGas_->RT();

                pmu[facei] = mixture_.CanteraTransport()->viscosity();

                palpha[facei] = mixture_.CanteraTransport()->thermalConductivity()/(CanteraGas_->cp_mass());

                if (mixture_.transportModelName() == "UnityLewis")
                {
                    forAll(rhoD_, i)
                    {
                        rhoD_[i].boundaryFieldRef()[patchi][facei] = palpha[facei];
                    }
                }
                else
                {
                    mixture_.CanteraTransport()->getMixDiffCoeffsMass(dTemp_.begin());

                    CanteraGas_->getEnthalpy_RT(hrtTemp_.begin());
                    const scalar RT = constant::physicoChemical::R.value()*1e3*pT[facei];
                    forAll(rhoD_, i)
                    {
                        rhoD_[i].boundaryFieldRef()[patchi][facei] = dTemp_[i];

                        hai_[i].boundaryFieldRef()[patchi][facei] = hrtTemp_[i]*RT/CanteraGas_->molecularWeight(i);
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
