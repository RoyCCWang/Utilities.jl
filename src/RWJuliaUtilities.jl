
module RWJuliaUtilities

using FFTW
import Random
using LinearAlgebra

import Printf

#import HCubature


import GSL
import SpecialFunctions

import Interpolations
#import StaticArrays
import StatsFuns

import Test

include("declarations.jl")

include("matrix.jl")
include("approximations.jl")
include("geometric_transforms.jl")
include("signals.jl")
include("boundary.jl")
include("patterns.jl")
include("bounds.jl")
include("combinatorics.jl")
include("template_signals.jl")
include("random_matrices.jl")
include("data_structures.jl")
include("simplex.jl")
include("function_interpolation.jl")
include("numerical.jl")


export  isnumericallyclose, # bounds.jl
        isnumericallyin,
        makenotclose!,

        # combinatorics
        coprime,
        EuclideanAlgotihm,

        # patterns.jl
        extractselectcols,
        finddiscrepancy,
        forcesymmetric,
        forcenegativesymmetric,
        coordinatetolattice1D,
        lattice1Dtocoordinate,
        ranges2collection,

        # boundary.jl
        extractframecenters,
        padzeros,
        evalzeropaddedsignal,


        # signals.jl
        rad2freq,
        freq2rad,
        computeDTFTviaformula,
        getdiracdeltaseq,
        spectralinversion,
        upsample,
        downsample,
        linearconv,
        forcesamesizeafterlinearconv,
        fastlinearconv,
        circconv,

        # random_matrices.jl
        generaterandomsym,
        generaterandomskew,
        generaterandomortho,
        generateStiefel,
        generaterandomhankel,
        generateindefmatrix,
        generateposdefmatrix,

        # template_signals.jl
        rectfunc,
        raisedcosinefunction,
        cosinetransitionlowpassfunction,
        cosinetransitionhighpassfunction,
        cosinetransitionbandpassfunction,

        # geometric_transforms.jl
        shiftthenstretch,
        shiftthencompress,
        evalstretchedsignal,
        convertcompactdomain,

        # approximations.jl
        smoothmin,
        smoothmax,
        smoothextrema,


        # matrix.jl
        forcesymmetric,
        naivesqrtpsdmatrix,
        array2matrix,
        makeblockdiagonalmatrix,
        invertblockmatrix,
        blockmatrixmatrixproduct,
        blockmatrixvectorproduct,

        # data_structures.jl
        packbcW,
        unpackbcW,
        parseCholeskyType,
        packCholeskyType,

        # simplex.jl
        simplexinteriortoEuclidean,
        Euclideantosimplex,

        # function_interpolation
        setupcubicitp,

        # numerical
        logsumexp
end
