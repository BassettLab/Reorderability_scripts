## Packages

using Pkg
using Eirene
using MAT
using DataFrames
using JLD
using Plots
using Measures
using StatsPlots
using RCall
using Statistics
using StatsBase
using HypothesisTests
using LinearAlgebra
using MultivariateStats

# Call R functions
R"library(TDA)"
R"source('~/Dropbox/Top Sim and Homog/Scripts/bottleneck_computations_functions2.R')"
R"library(igraph)"
R"library(ggplot2)"
R"source('~/Dropbox/Top Sim and Homog/Scripts/local_network_functions.R')"





## Functions

function toWeightedAdj_byweight(network_und,node_weights)
    #Assume that HIGHER node weight means born earlier-- this means we
    # want to keep the MINIMUM
    nNodes = length(node_weights)
    t0 = network_und.*node_weights
    t1 = deepcopy(t0)
    for i in collect(1:nNodes)
        for j in collect(i:nNodes)
            t1[i,j] = minimum([t0[i,j],t0[j,i]])
            t1[j,i] = t1[i,j]
        end
    end
    t1
end

# To calculate the bettiCurves from the barcodes
function bettiCurveFromBarcode(barcode_array,nNodes,nmats,maxDim)

    nNodes = Int(nNodes)
    nmats = Int(nmats)
    maxDim = Int(maxDim)
    bettiBar = zeros(nmats,maxDim)
    bettiCurve = zeros(nmats,nNodes+1,maxDim)
    birthCurve = zeros(nmats,nNodes,maxDim)
    deathCurve = zeros(nmats,nNodes,maxDim)


    for dimn in collect(1:maxDim)
        dimn = Int(dimn)

        for matn in collect(1:nmats)
            matn = Int(matn)
            bb = 0
            currentCurve = barcode_array[matn,:]
            currentCurveDim = currentCurve[dimn]
            for barn in collect(1:size(currentCurveDim,1))


                # Add to birth curve
                birthCurve[matn,Int(currentCurveDim[barn,1]),dimn] = birthCurve[matn,Int(currentCurveDim[barn,1]),dimn] .+1


                if currentCurveDim[barn,2]>nNodes

                    bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(nNodes+1),dimn] = bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(nNodes+1),dimn] .+1
                    bb = bb+(nNodes+1-currentCurveDim[barn,1])
                    else

                    bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(currentCurveDim[barn,2]),dimn] = bettiCurve[matn,Int(currentCurveDim[barn,1]):Int(currentCurveDim[barn,2]),dimn].+1
                    deathCurve[matn,Int(currentCurveDim[barn,2]),dimn] = deathCurve[matn,Int(currentCurveDim[barn,2]),dimn] .+1
                    bb = bb+(currentCurveDim[barn,2] - currentCurveDim[barn,1])

                end
            end

            bettiBar[matn,dimn] = deepcopy(bb)
        end
    end

    return bettiCurve, birthCurve, deathCurve, bettiBar
end


## Calculating edge density
function calculateEdgeQuantity(jadj_array)
    nmats = size(jadj_array,3)
    nNodes = size(jadj_array,1)
    edgeQuantity = zeros(nmats,nNodes)

    # Loop through to compute number of edges at each filtration step
    for matn in collect(1:nmats)
        adj = jadj_array[:,:,matn];
        for noden in collect(1:nNodes)
            edgesAdded = (LinearIndices(adj))[findall(adj.<= noden)]
            edgeQuantity[Int(matn),Int(noden)] = length(edgesAdded)/2
        end
    end
    edgeQuantity = edgeQuantity.-35   # get rid of self-loops = 0

end

function calculateDegreesFiltration(jadj_array,nNodes)
    nmats = size(jadj_array,3)
    nNodes = size(jadj_array,1)
    degree_array = zeros(nNodes,nNodes,nmats)

    jadj_array[jadj_array.==(nNodes*2)] .= 0


    # Loop through to compute number of edges at each filtration step
    for matn in collect(1:nmats)
        adj = jadj_array[:,:,matn];
        for noden in collect(1:nNodes)
            adj_i = deepcopy(adj)
            adj_i[adj_i.>noden] .= 0
            adj_i[adj_i.>0] .= 1
            degree_array[:,Int(noden),Int(matn)] = sum(adj_i,dims = 1)
        end
    end

    return degree_array

end


## Making weights from order
function orderToWeights(s_0_array,nNodes)
        s_wei_array = zeros(size(s_0_array))
    for i0 in collect(1:size(s_0_array,1))

        for n0 in collect(1:nNodes)
            n0 = Int(n0)
            s_wei_array[i0,n0] = indexin([n0],s_0_array[i0,:])[1]
        end
    end
    s_wei_array
end

function calculateTheoreticaMaxDistance(s_wei_array,n)
    theo_dist= Vector{Float64}(undef,n)   #10100 or 24160
    for i1 in collect(1:size(s_wei_array,1))
        norm1 = norm(s_wei_array[1,:].-s_wei_array[i1,:],Inf)
        theo_dist[i1] = norm1
    end
    return theo_dist
end

## Making Diagrams from barcode for TDA in R
function makeDiagramFromBarcode(barcode_array,graphn,maxDim)

    diag = []
    nNodes = 70


    for d in collect(1:maxDim)

        bd = barcode_array[graphn,d]
        infs = findall(bd[:,2] .> nNodes)
        bd[infs,2] .= nNodes+1
        nbars = size(bd,1)

        if nbars>0
            bda = hcat(d*ones(nbars),bd)
        else
            bda = [d 0 0]
        end




        if isequal(diag,[])
            diag = bda
            else
            diag = vcat(diag,bda)
        end

    end
    diag
end




function computeBNDistances_glob(barcode_array)

    nReps = 101
    nGraphs = 100
    maxDim = 3

    graphOriginals = collect(1:101:10100)
    distanceBN_array = zeros(10100,3)


    for nG in graphOriginals

        diagO = makeDiagramFromBarcode(barcode_array,nG,maxDim)

        for nR in collect(1:(nReps-1))
            diagR = makeDiagramFromBarcode(barcode_array,(nG+nR),maxDim)

            for nD in collect(1:3)
                distanceBN_dimn = R"computeBNDistance($diagO,$diagR,$nD)"
                distanceBN_array[(nG+nR),nD] = distanceBN_dimn
            end
        end
    end

    return distanceBN_array
end



function computeBNDistances_globSampled(barcode_array)

    nReps = 101
    nGraphs = 100
    maxDim = 3
    runs = collect(1:100)

    graphOriginals = collect(1:101:10100)
    distanceBN_array = zeros(10000,3)


    for nG in collect(1:nGraphs)

        a = (nG-1)*nReps + 1

        for runi in runs

        p = sample(collect(a:(a+100)),2,replace = false)
        diaga = makeDiagramFromBarcode(barcode_array,p[1],maxDim)
        diagb = makeDiagramFromBarcode(barcode_array,p[2],maxDim)

            for nD in collect(1:3)
                distanceBN_dimn = R"computeBNDistance($diaga,$diagb,$nD)"
                distanceBN_array[((nG-1)*(nReps-1) + runi),nD] = distanceBN_dimn
            end
        end

    end

    return distanceBN_array
end




function computeBNDistances_Sampled(barcode_array,nGraphs)


    maxDim = 3
    distanceBN_array = zeros(10000,3)

    runs = collect(1:size(distanceBN_array)[1])

    for runi in runs
        p = sample(collect(1:nGraphs),2,replace = false)
        diaga = makeDiagramFromBarcode(barcode_array,p[1],maxDim)
        diagb = makeDiagramFromBarcode(barcode_array,p[2],maxDim)

        for nD in collect(1:3)
            distanceBN_dimn = R"computeBNDistance($diaga,$diagb,$nD)"
            distanceBN_array[runi,nD] = distanceBN_dimn
        end

    end

    return distanceBN_array
end


function computeBNDistances_local(barcode_array)

    nReps = 2416
    nGraphs = 10
    maxDim = 3
    nNodes = 70

    graphOriginals = collect(1:2416:24160)
    distanceBN_array = zeros(24160,3)

    for nG in graphOriginals
        diagO = makeDiagramFromBarcode(barcode_array,nG,maxDim)
        println(nG)
        for nR in collect(1:(nReps-1))
            diagR = makeDiagramFromBarcode(barcode_array,(nG+nR),maxDim)

            for nD in collect(1:3)
                distanceBN_dimn = R"computeBNDistance($diagO,$diagR,$nD)"
                distanceBN_array[(nG+nR),nD] = distanceBN_dimn
            end
        end
    end

    return distanceBN_array
end


function plotBarcode(allPIs,nNodes,graphN,maxDim,fontSize)

    nNodes = Int(nNodes)
    graphn = Int(graphN)
    maxDim = Int(maxDim)
    counter1 = 0
    pbar = plot(1:6,zeros(6),c=:black)

    colors = [:blue :green :red]
    for dim in collect(1:maxDim)

        barn = barcode_array[graphN, dim]
        barn = barn[sortperm(barn[:,1]),:]

        nbars = size(barn)[1]


        for cntr1 in collect(1:nbars)
            birth = barn[cntr1,1]
            death = barn[cntr1,2]

            plot!([birth, death],[cntr1+counter1, cntr1+counter1],c=colors[dim], legend = false,
                            xlim = (0,nNodes), ytickfont = font(fontSize), xtickfont = font(fontSize))
        end

        display(pbar)



        counter1 = counter1+nbars
    end

    return pbar
end


function strictLTvector(x)
    xVec = [x[i,j] for j in 1:size(x)[2]-1 for i in j+1:size(x)[1]]
    return xVec
end

function strictLTvector3D(x)
    xVec = [x[i,j,d] for d in 1:size(x)[3] for j in 1:size(x)[2]-1 for i in j+1:size(x)[1]]
    return xVec
end
