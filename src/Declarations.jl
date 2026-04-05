


function HomalgRingOfIntegersInSingular end; export HomalgRingOfIntegersInSingular

function HomalgFieldOfRationalsInSingular end; export HomalgFieldOfRationalsInSingular

function AssignGeneratingVariables end; export AssignGeneratingVariables

function HasHasInvariantBasisProperty end; export HasHasInvariantBasisProperty

HasHasInvariantBasisProperty(::TypeOfRingForHomalg) = false
HasHasInvariantBasisProperty(::Union{typeof(QQ), typeof(ZZ)}) = true

function HasInvariantBasisProperty end; export HasInvariantBasisProperty

HasInvariantBasisProperty(::Union{typeof(QQ), typeof(ZZ)}) = true;

function RandomMatrix end; export RandomMatrix

function Indeterminates end; export Indeterminates
