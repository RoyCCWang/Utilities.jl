
module Utilities

using FFTW
import Random
using LinearAlgebra
import PyPlot
import Printf

import HCubature

import Distributions
import GSL
import SpecialFunctions

import Interpolations
#import StaticArrays
import StatsFuns

import Test

include("declarations.jl")
include("visualize_DSP.jl")
include("visualize.jl")
include("visualize_2D.jl")
include("matrix.jl")
include("approximations.jl")
include("geometric_transforms.jl")
include("signals.jl")
include("boundary.jl")
include("patterns.jl")
include("bounds.jl")
include("combinatorics.jl")
include("template_signals.jl")
include("drawRV.jl")
include("probability.jl")
include("probability_metrics.jl")
include("data_structures.jl")
include("simplex.jl")
include("function_interpolation.jl")
include("numerical.jl")
include("synthetic_densities.jl")

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

        # synthetic_densities.jl
        generaterandomGMM,
        getGMMdist,

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


        # # visualize.jl
        visualizsignals,
        visualizefbmagrsp,
        plotmagnitudersp,
        plotphasersp,
        getfreqrsp

        # drawRV.jl
        drawfromuniformindexinterval,
        generaterandomsym,
        generaterandomskew,
        generaterandomortho,
        generateStiefel,
        generaterandomhankel,
        generateindefmatrix,
        generateposdefmatrix,
        probit,
        Gaussiancopula,
        evalGaussiancopula,
        drawstdnormalcdf,
        drawstdnormalcdf!

        # # visualize_2D.jl
        visualizemeshgridpcolor,
        plot2Dhistogram,

        # probability.jl
        getnormalizedensity,
        evalstdnormalcdf,
        evalstdnormalpdf,
        getMVNmarginalparams,

        # probability_metrics.jl
        evalKLintegral,

        # data_structures.jl
        packbcW,
        unpackbcW,
        parseCholeskyType,
        packCholeskyType,

        # simplex.jl
        simplexinteriortoEuclidean,
        Euclideantosimplex,

        # # visualize.jl
        plothistogram,

        # function_interpolation
        setupcubicitp,

        # numerical
        logsumexp
end
