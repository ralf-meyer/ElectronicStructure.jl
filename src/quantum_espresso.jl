using ASEconvert
using PythonCall
using MPI
using LinearAlgebra
using Unitful
using UnitfulAtomic


Base.@kwdef struct QeCalculator <: AbstractCalculator
    # TODO QE setup parameters
end

Base.@kwdef struct QeParameters <: AbstractParameters
    # TODO Keywords currently based on ASE
    system::AbstractSystem
    ecutwfc     = 40
    conv_thr    = 1e-11
    tstress     = true
    tprnfor     = true
    smearing    = "gaussian"
    mixing_mode = "plain"
    mixing_beta = 0.7
    mixing_ndim = 10
    kpts        = (1, 1, 1)
    occupations = "smearing"
    degauss     = 0.01
    input_dft   = "pbe"
    electron_maxstep  = 100
    pseudopotentials  = Dict{String,String}()
    extra_parameter   = Dict{Symbol,Any}()
    working_directory = mktempdir(pwd())
    n_mpi_procs       = MPI.Comm_size(MPI.COMM_WORLD)
    n_threads         = BLAS.get_num_threads()
end



function convert(::Type{QeParameters}, params::DftkParameters)
    # Convert DFTK parameters to QE parameters
    #
    # keep in mind the unit conversion
    error("TODO")
end


struct QeState <: AbstractState
    params::QeParameters
    ase_atoms::Py
end

function QeState(params::QeParameters)
    ase_atoms = convert_ase(params.system)
    ase_atoms.calc = pyimport("ase.calculators.espresso").Espresso(;
        label="espresso",
        ecutwfc=params.ecutwfc,
        conv_thr=params.conv_thr,
        tstress=params.tstress,
        tprnfor=params.tprnfor,
        smearing=params.smearing,
        mixing_mode=params.mixing_mode,
        mixing_beta=params.mixing_beta,
        mixing_ndim=params.mixing_ndim,
        kpts=params.kpts,
        occupations=params.occupations,
        degauss=params.degauss,
        input_dft=params.input_dft,
        electron_maxstep=params.electron_maxstep,
        pseudopotentials=params.pseudopotentials,
        params.extra_parameter...
    )
    QeState(params, ase_atoms)
end

function calculate(calc::QeCalculator, params::QeParameters)
    calculate(calc, QeState(params))
end

function calculate(::QeCalculator, state::QeState)
    n_mpi_procs = state.params.n_mpi_procs
    qe_path = ENV["ESPRESSO_BIN"]
    qe_command = "mpirun -np $n_mpi_procs $qe_path/pw.x -in PREFIX.pwi > PREFIX.pwo"
    state.ase_atoms.calc.command = qe_command
    withenv("OMP_NUM_THREADS" => state.params.n_threads) do
        state.ase_atoms.get_potential_energy()
    end
end

function energy(state::QeState)
    austrip(pyconvert(AbstractFloat, state.ase_atoms.get_potential_energy())u"eV")
end
