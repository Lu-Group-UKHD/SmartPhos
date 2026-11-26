# Get Ratio Matrix of Phosphoproteome Data

`getRatioMatrix` calculates the ratio matrix of phosphoproteome data
from a `MultiAssayExperiment` object.

## Usage

``` r
getRatioMatrix(maeData, normalization = FALSE, getAdjustedPP = FALSE)
```

## Arguments

- maeData:

  A `MultiAssayExperiment` object containing phosphoproteome and full
  proteome data.

- normalization:

  A `logical` value indicating whether to perform normalization. Default
  is `FALSE`.

- getAdjustedPP:

  A `logical` value indicating whether to use adjusted phosphoproteome
  data. Default is `FALSE`.

## Value

A numeric `matrix` representing the ratio of intensity of PP
(phosphoproteome) data to FP (full proteome) data.

## Examples

``` r
# Load multiAssayExperiment object
data("dia_example")
# Call the function
getRatioMatrix(dia_example, normalization = TRUE)
#>      Phospho_1stCrtl_0min_rep2 Phospho_1stCrtl_0min_rep3
#> s13                         NA                        NA
#> s16                 -3.7643939                        NA
#> s36                         NA                        NA
#> s37                         NA                        NA
#> s40                         NA                        NA
#> s57                 -4.3233417                -3.3027041
#> s58                 -4.3233417                -3.3027041
#> s59                 -4.3233417                -3.3027041
#> s60                 -4.9211318                -1.8669683
#> s64                         NA                        NA
#> s65                         NA                -3.1975354
#> s79                         NA                        NA
#> s83                         NA                        NA
#> s93                         NA                        NA
#> s94                  2.2532894                        NA
#> s114                        NA                        NA
#> s128                        NA                        NA
#> s129                        NA                        NA
#> s136                        NA                        NA
#> s142                        NA                        NA
#> s144                        NA                        NA
#> s145                        NA                 8.8920740
#> s147                        NA                        NA
#> s150                -1.5001613                        NA
#> s151                        NA                        NA
#> s153                -1.7884873                        NA
#> s154                        NA                        NA
#> s155                -0.2206020                        NA
#> s158                        NA                        NA
#> s162                        NA                        NA
#> s167                        NA                -1.6745031
#> s173                -1.9308668                 1.6516496
#> s174                -4.3686259                        NA
#> s176                        NA                        NA
#> s177                        NA                        NA
#> s179                        NA                        NA
#> s180                        NA                        NA
#> s181                        NA                        NA
#> s188                        NA                -1.2480599
#> s197                        NA                        NA
#> s198                        NA                        NA
#> s199                        NA                        NA
#> s212                        NA                 0.1952259
#> s219                        NA                        NA
#> s220                        NA                        NA
#> s221                        NA                        NA
#> s234                        NA                        NA
#> s237                        NA                        NA
#> s240                        NA                        NA
#> s244                        NA                        NA
#> s246                        NA                        NA
#> s265                 3.5204058                        NA
#> s270                -8.2796228                -7.0913119
#> s271                -8.2796228                -7.0913119
#> s274                -6.9236588                -8.3827098
#> s278                        NA                        NA
#> s283                        NA                -6.3197134
#> s287                        NA                        NA
#> s299                        NA                -3.2379105
#> s302                        NA                        NA
#> s304                        NA                        NA
#> s307                        NA                -3.0786511
#> s311                        NA               -10.0262252
#> s313                        NA                        NA
#> s322                        NA                        NA
#> s325                        NA                        NA
#> s330                        NA                        NA
#> s341                        NA                        NA
#> s347                        NA                        NA
#> s348                        NA                        NA
#> s356                -1.5070661                        NA
#> s362                 0.3978409                -0.1927442
#> s374                        NA                        NA
#> s376                -0.3729111                -1.1199006
#> s377                -7.2027484                        NA
#> s381                -3.7637242                -3.8044981
#> s382                        NA                        NA
#> s393                 0.2304416                -0.4156391
#> s395                        NA                        NA
#> s396                        NA                 3.6908447
#> s401                        NA                        NA
#> s416                        NA                        NA
#> s418                        NA                        NA
#> s430                        NA                        NA
#> s436                -7.0904517                        NA
#> s444                -2.1876346                        NA
#> s448                        NA                -9.7025047
#> s451                        NA                -4.0506228
#> s452                        NA                 1.0338960
#> s454                        NA                        NA
#> s457                        NA                        NA
#> s470                        NA                -2.7403064
#> s472                        NA                        NA
#> s474                 0.5698129                -0.0758101
#> s482                        NA                        NA
#> s488                        NA                        NA
#> s499                        NA                        NA
#> s500                 1.3884624                 0.7454984
#>      Phospho_1stCrtl_100min_rep2 Phospho_1stCrtl_100min_rep3
#> s13                           NA                          NA
#> s16                           NA                          NA
#> s36                           NA                          NA
#> s37                           NA                          NA
#> s40                  2.414032293                  3.59138064
#> s57                 -3.446662944                 -2.51928708
#> s58                 -3.446662944                 -2.51928708
#> s59                 -3.446662944                 -2.51928708
#> s60                 -3.446662944                 -1.45078302
#> s64                           NA                          NA
#> s65                           NA                          NA
#> s79                           NA                 -1.41597841
#> s83                           NA                          NA
#> s93                           NA                          NA
#> s94                           NA                          NA
#> s114                          NA                          NA
#> s128                          NA                 -9.78903466
#> s129                          NA                 -9.78903466
#> s136                          NA                          NA
#> s142                          NA                          NA
#> s144                          NA                          NA
#> s145                -4.533103313                 -1.12185247
#> s147                          NA                          NA
#> s150                          NA                          NA
#> s151                          NA                          NA
#> s153                          NA                          NA
#> s154                -2.489925626                          NA
#> s155                          NA                  0.48252290
#> s158                          NA                 -2.03910048
#> s162                          NA                          NA
#> s167                          NA                 -4.54359072
#> s173                -2.531618195                 -2.10220759
#> s174                          NA                          NA
#> s176                          NA                          NA
#> s177                          NA                          NA
#> s179                          NA                          NA
#> s180                          NA                          NA
#> s181                -4.343526443                          NA
#> s188                -1.095991743                          NA
#> s197                          NA                          NA
#> s198                 2.212898838                          NA
#> s199                          NA                          NA
#> s212                 1.483997733                  0.26346771
#> s219                          NA                          NA
#> s220                          NA                          NA
#> s221                          NA                  0.09330300
#> s234                          NA                  9.54959755
#> s237                          NA                          NA
#> s240                 3.882016582                          NA
#> s244                          NA                 -2.94990304
#> s246                          NA                          NA
#> s265                          NA                          NA
#> s270               -10.341975085                 -9.79108232
#> s271               -10.341975085                 -9.79108232
#> s274                -8.936996355                 -8.43425891
#> s278                          NA                          NA
#> s283                          NA                          NA
#> s287                          NA                          NA
#> s299                -3.730188485                 -5.25817801
#> s302                          NA                          NA
#> s304                          NA                          NA
#> s307                -4.430467039                 -2.22740000
#> s311                          NA                          NA
#> s313                          NA                          NA
#> s322                 3.157243008                          NA
#> s325                          NA                          NA
#> s330                          NA                          NA
#> s341                          NA                          NA
#> s347                          NA                          NA
#> s348                          NA                          NA
#> s356                          NA                          NA
#> s362                 0.004915533                  0.07296485
#> s374                          NA                          NA
#> s376                -0.943445996                 -1.47253029
#> s377                          NA                          NA
#> s381                -2.567132319                 -3.24466915
#> s382                          NA                          NA
#> s393                 0.001857241                 -0.15719872
#> s395                          NA                          NA
#> s396                          NA                  3.14368385
#> s401                          NA                          NA
#> s416                          NA                          NA
#> s418                          NA                          NA
#> s430                          NA                          NA
#> s436                -7.254007291                          NA
#> s444                          NA                          NA
#> s448                          NA                          NA
#> s451                          NA                          NA
#> s452                          NA                 -0.74753146
#> s454                          NA                          NA
#> s457                          NA                          NA
#> s470                 0.010617871                          NA
#> s472                          NA                          NA
#> s474                 1.716430997                  0.46489341
#> s482                -3.728395761                 -0.36131445
#> s488                          NA                          NA
#> s499                          NA                  2.40259179
#> s500                 0.748722186                  0.83472183
#>      Phospho_EGF_10min_rep2 Phospho_EGF_10min_rep3 Phospho_EGF_100min_rep2
#> s13            -7.025506174                     NA                      NA
#> s16            -5.556174953                     NA              -3.8222037
#> s36                      NA                     NA                      NA
#> s37                      NA              0.8709432                      NA
#> s40                      NA                     NA               2.5930917
#> s57            -2.634147566             -4.5640618              -2.2585646
#> s58            -2.634147566             -4.5640618              -2.2585646
#> s59            -2.634147566             -4.5640618              -2.2585646
#> s60            -2.634147566             -5.2063405              -2.2585646
#> s64             1.631715360                     NA                      NA
#> s65            -0.006010127                     NA                      NA
#> s79                      NA                     NA                      NA
#> s83             4.083505008                     NA                      NA
#> s93                      NA                     NA                      NA
#> s94                      NA                     NA                      NA
#> s114           -5.960754232                     NA                      NA
#> s128                     NA                     NA                      NA
#> s129                     NA                     NA                      NA
#> s136                     NA                     NA                      NA
#> s142           -6.895033509                     NA                      NA
#> s144                     NA                     NA                      NA
#> s145                     NA                     NA              -2.2053733
#> s147                     NA              4.7138385                      NA
#> s150           -0.088360496                     NA              -1.7554356
#> s151                     NA                     NA              -1.7371280
#> s153                     NA                     NA                      NA
#> s154           -0.030058853                     NA                      NA
#> s155                     NA                     NA              -0.8632875
#> s158                     NA                     NA              -4.4545300
#> s162            1.558051503              4.9800398                      NA
#> s167                     NA                     NA                      NA
#> s173           -2.715464977              0.9857779              -2.9365945
#> s174                     NA                     NA                      NA
#> s176                     NA             -6.4248242             -10.1680540
#> s177                     NA             -6.4248242              -8.7199236
#> s179                     NA             -1.3044115                      NA
#> s180                     NA                     NA                      NA
#> s181                     NA                     NA                      NA
#> s188           -1.088815486             -0.6549753                      NA
#> s197           -5.282908049                     NA                      NA
#> s198                     NA                     NA                      NA
#> s199                     NA                     NA                      NA
#> s212            0.264601131              0.3043399               0.4239538
#> s219            0.653740249                     NA                      NA
#> s220            4.296131588                     NA                      NA
#> s221                     NA                     NA               0.5698976
#> s234           10.101989461                     NA               0.2201499
#> s237                     NA                     NA                      NA
#> s240                     NA                     NA                      NA
#> s244                     NA                     NA                      NA
#> s246                     NA                     NA                      NA
#> s265                     NA                     NA                      NA
#> s270                     NA            -10.5271969                      NA
#> s271                     NA            -10.5271969                      NA
#> s274           -8.275496761             -7.8471385              -7.6725834
#> s278           -2.264472152                     NA              -3.7080424
#> s283                     NA                     NA                      NA
#> s287           -0.737231627                     NA              -1.2037960
#> s299           -1.745238954             -3.6697727              -5.2699696
#> s302                     NA             -0.9554264                      NA
#> s304                     NA                     NA                      NA
#> s307           -3.248705393             -3.3531115              -3.8507753
#> s311                     NA                     NA                      NA
#> s313                     NA                     NA                      NA
#> s322            1.644023592             -1.9626936                      NA
#> s325                     NA                     NA               0.8947103
#> s330                     NA                     NA                      NA
#> s341                     NA                     NA              -5.0873552
#> s347           -0.637347404                     NA                      NA
#> s348                     NA                     NA                      NA
#> s356                     NA                     NA                      NA
#> s362           -0.837153924                     NA              -0.1999119
#> s374                     NA                     NA                      NA
#> s376           -0.986094449             -2.0678469              -1.2613732
#> s377                     NA                     NA              -0.7666310
#> s381           -4.407529659             -4.3554403              -5.0762369
#> s382                     NA                     NA                      NA
#> s393           -0.414638314             -0.2392067              -0.2934095
#> s395                     NA                     NA                      NA
#> s396            2.835465953                     NA                      NA
#> s401                     NA                     NA                      NA
#> s416                     NA                     NA                      NA
#> s418                     NA                     NA              -0.2748643
#> s430           -2.220519607                     NA                      NA
#> s436           -7.509013050             -7.1500188                      NA
#> s444                     NA                     NA                      NA
#> s448                     NA                     NA                      NA
#> s451           -6.090285030             -6.4500009                      NA
#> s452                     NA                     NA                      NA
#> s454                     NA                     NA                      NA
#> s457                     NA                     NA                      NA
#> s470                     NA             -3.4419921                      NA
#> s472                     NA                     NA                      NA
#> s474            0.095250177              0.3781762               0.1442772
#> s482                     NA             -1.6690894              -1.2880232
#> s488                     NA                     NA              -7.9628736
#> s499                     NA                     NA                      NA
#> s500            0.149844819              0.7993372               0.3023715
#>      Phospho_EGF_100min_rep3 Phospho_HGF_10min_rep2 Phospho_HGF_10min_rep3
#> s13                       NA                     NA                     NA
#> s16                       NA            -3.94015233                     NA
#> s36                       NA             2.15897383                     NA
#> s37                       NA                     NA                     NA
#> s40                       NA             2.44515585                     NA
#> s57              -4.18897750            -4.10933733             -2.9207293
#> s58              -4.18897750            -4.10933733             -2.9207293
#> s59              -4.18897750            -4.10933733             -2.9207293
#> s60              -4.18897750            -4.91682017             -2.9207293
#> s64                       NA                     NA              1.0668834
#> s65                       NA                     NA                     NA
#> s79                       NA                     NA                     NA
#> s83                       NA                     NA                     NA
#> s93              -2.28322345                     NA                     NA
#> s94                       NA                     NA                     NA
#> s114                      NA                     NA                     NA
#> s128                      NA                     NA                     NA
#> s129                      NA                     NA                     NA
#> s136             14.31235974                     NA                     NA
#> s142                      NA            -5.41535671                     NA
#> s144                      NA                     NA                     NA
#> s145                      NA                     NA                     NA
#> s147                      NA                     NA                     NA
#> s150                      NA                     NA                     NA
#> s151             -3.14211361                     NA                     NA
#> s153             -1.78585633                     NA                     NA
#> s154                      NA             1.29184206                     NA
#> s155                      NA                     NA                     NA
#> s158             -3.50322539                     NA                     NA
#> s162                      NA                     NA                     NA
#> s167             -3.99095421             6.62157929                     NA
#> s173             -1.35581798            -1.11159628              1.7085950
#> s174                      NA                     NA                     NA
#> s176                      NA                     NA                     NA
#> s177                      NA                     NA                     NA
#> s179                      NA                     NA                     NA
#> s180                      NA                     NA                     NA
#> s181                      NA                     NA                     NA
#> s188             -0.92222013                     NA             -0.4928061
#> s197                      NA            -5.05887765                     NA
#> s198                      NA                     NA                     NA
#> s199                      NA            -1.09900797                     NA
#> s212                      NA             0.57234627              0.6497259
#> s219                      NA                     NA                     NA
#> s220                      NA                     NA                     NA
#> s221                      NA                     NA                     NA
#> s234                      NA                     NA                     NA
#> s237                      NA                     NA              1.7622516
#> s240                      NA                     NA                     NA
#> s244                      NA            -3.26628311                     NA
#> s246                      NA                     NA                     NA
#> s265              0.07735383                     NA                     NA
#> s270            -10.14583417            -7.27460097             -9.8822933
#> s271            -10.14583417            -7.27460097             -9.8822933
#> s274             -7.67927731            -8.12701776             -7.6784462
#> s278                      NA                     NA                     NA
#> s283                      NA                     NA                     NA
#> s287                      NA                     NA                     NA
#> s299             -4.76153981            -5.31066046             -2.5101768
#> s302                      NA                     NA                     NA
#> s304              1.11599956                     NA                     NA
#> s307                      NA                     NA                     NA
#> s311                      NA                     NA                     NA
#> s313                      NA                     NA                     NA
#> s322             -2.51800318            -1.99753838                     NA
#> s325                      NA                     NA                     NA
#> s330                      NA             3.55241926              2.2394399
#> s341                      NA                     NA                     NA
#> s347                      NA                     NA                     NA
#> s348                      NA             3.24974688                     NA
#> s356                      NA                     NA                     NA
#> s362             -0.18130222            -0.08825152              0.4569485
#> s374                      NA                     NA                     NA
#> s376             -2.34717624             1.89361401             -0.4842351
#> s377                      NA                     NA             -3.0284911
#> s381             -4.63558862            -4.86594529             -4.2492993
#> s382                      NA                     NA                     NA
#> s393             -0.46275198             0.21203169             -0.1676951
#> s395             -3.69051651                     NA                     NA
#> s396              2.94197556             2.87511952              2.7616507
#> s401                      NA                     NA                     NA
#> s416                      NA                     NA                     NA
#> s418                      NA                     NA                     NA
#> s430                      NA                     NA                     NA
#> s436                      NA            -7.54740416                     NA
#> s444                      NA                     NA                     NA
#> s448                      NA                     NA                     NA
#> s451                      NA                     NA                     NA
#> s452                      NA                     NA                     NA
#> s454              0.64227255                     NA                     NA
#> s457                      NA                     NA              0.5700706
#> s470                      NA                     NA             -0.9069880
#> s472              6.47432216                     NA                     NA
#> s474             -0.25805015             1.05040336              0.3613268
#> s482                      NA                     NA                     NA
#> s488                      NA                     NA                     NA
#> s499                      NA                     NA                     NA
#> s500             -1.06057119             0.56263283                     NA
#>      Phospho_HGF_100min_rep2 Phospho_HGF_100min_rep3 Phospho_1stCrtl_0min_rep1
#> s13                       NA              -6.0252423                        NA
#> s16                       NA                      NA                        NA
#> s36                       NA                      NA                        NA
#> s37                       NA                      NA                        NA
#> s40                       NA                      NA                        NA
#> s57              -1.68874750              -1.1364318               -2.68632848
#> s58              -1.68874750              -1.1364318               -2.68632848
#> s59              -1.68874750              -1.1364318               -2.68632848
#> s60              -2.48774011              -1.1364318               -2.68632848
#> s64                       NA                      NA                        NA
#> s65                       NA                      NA                        NA
#> s79                       NA                      NA                        NA
#> s83                       NA                      NA                        NA
#> s93                       NA                      NA                        NA
#> s94                       NA                      NA                        NA
#> s114                      NA                      NA                        NA
#> s128                      NA                      NA                        NA
#> s129                      NA                      NA                        NA
#> s136                      NA                      NA                        NA
#> s142                      NA                      NA                        NA
#> s144              2.36849245                      NA                        NA
#> s145                      NA                      NA                        NA
#> s147                      NA               0.4245685                        NA
#> s150                      NA                      NA                        NA
#> s151                      NA                      NA               -1.77563524
#> s153                      NA                      NA                        NA
#> s154                      NA                      NA                        NA
#> s155                      NA                      NA                        NA
#> s158                      NA                      NA                        NA
#> s162                      NA               3.2412030                        NA
#> s167                      NA                      NA                        NA
#> s173             -0.73211692              -2.3106651               -1.67778572
#> s174                      NA                      NA                        NA
#> s176             -5.79107561                      NA                        NA
#> s177             -5.79107561                      NA                        NA
#> s179                      NA                      NA                        NA
#> s180                      NA              -3.0472205                        NA
#> s181                      NA                      NA                        NA
#> s188             -0.83201254              -0.3813293                0.82120229
#> s197                      NA                      NA                        NA
#> s198                      NA                      NA                        NA
#> s199                      NA                      NA                        NA
#> s212              0.63430078               0.9853540                        NA
#> s219              2.08969581                      NA                        NA
#> s220                      NA                      NA                        NA
#> s221              1.04459780                      NA                        NA
#> s234                      NA                      NA               11.37116345
#> s237                      NA                      NA                3.11567330
#> s240                      NA                      NA                        NA
#> s244                      NA                      NA                        NA
#> s246                      NA                      NA                        NA
#> s265                      NA                      NA                        NA
#> s270                      NA                      NA                        NA
#> s271                      NA                      NA                        NA
#> s274             -7.59422976                      NA                        NA
#> s278             -4.10492801                      NA                        NA
#> s283                      NA                      NA                        NA
#> s287                      NA                      NA                        NA
#> s299             -5.72084076                      NA               -2.76913358
#> s302                      NA                      NA                        NA
#> s304                      NA                      NA                        NA
#> s307                      NA              -3.9611183                2.12882906
#> s311                      NA                      NA                        NA
#> s313                      NA              -4.3989730                        NA
#> s322                      NA                      NA                        NA
#> s325                      NA                      NA                        NA
#> s330                      NA                      NA                        NA
#> s341                      NA                      NA                        NA
#> s347                      NA                      NA                        NA
#> s348                      NA                      NA                        NA
#> s356                      NA                      NA                        NA
#> s362              0.34632966               0.6241446                0.09162791
#> s374             13.81234052                      NA                        NA
#> s376             -1.19834382              -0.7339090                        NA
#> s377                      NA                      NA                        NA
#> s381             -4.57087062              -3.7681520               -0.39839266
#> s382             -1.96234183                      NA                        NA
#> s393             -0.03704081              -0.1992859               -0.98039387
#> s395             -2.61268104                      NA                        NA
#> s396                      NA                      NA                        NA
#> s401                      NA               5.0001281                        NA
#> s416                      NA                      NA                        NA
#> s418              0.95099717                      NA                        NA
#> s430                      NA                      NA                        NA
#> s436             -7.51870156                      NA               -5.39444488
#> s444                      NA                      NA                        NA
#> s448                      NA                      NA               -6.68121086
#> s451                      NA                      NA                        NA
#> s452                      NA               0.1777189                        NA
#> s454                      NA                      NA                        NA
#> s457                      NA                      NA                        NA
#> s470             -0.30910237                      NA                        NA
#> s472                      NA                      NA                        NA
#> s474              0.25704761               0.6816826               -0.54158350
#> s482                      NA                      NA                        NA
#> s488             -8.21496276                      NA                        NA
#> s499                      NA                      NA                        NA
#> s500              1.25784743               1.2933145                1.51811754
#>      Phospho_1stCrtl_100min_rep1 Phospho_EGF_10min_rep1 Phospho_EGF_24h_rep1
#> s13                           NA                     NA         -10.07403491
#> s16                           NA                     NA                   NA
#> s36                           NA                     NA                   NA
#> s37                           NA                     NA          10.80789513
#> s40                           NA                     NA                   NA
#> s57                   3.19752774            -1.28417467          -0.10970776
#> s58                   3.19752774            -1.28417467          -0.10970776
#> s59                   3.19752774            -1.28417467          -0.10970776
#> s60                   3.19752774            -1.28417467          -0.10970776
#> s64                           NA                     NA                   NA
#> s65                           NA                     NA                   NA
#> s79                           NA                     NA                   NA
#> s83                           NA                     NA                   NA
#> s93                           NA                     NA                   NA
#> s94                           NA                     NA                   NA
#> s114                          NA                     NA                   NA
#> s128                          NA                     NA                   NA
#> s129                          NA                     NA                   NA
#> s136                          NA                     NA                   NA
#> s142                          NA                     NA                   NA
#> s144                          NA                     NA                   NA
#> s145                          NA                     NA                   NA
#> s147                          NA                     NA                   NA
#> s150                          NA                     NA                   NA
#> s151                          NA                     NA                   NA
#> s153                          NA                     NA                   NA
#> s154                          NA                     NA                   NA
#> s155                          NA                     NA                   NA
#> s158                          NA                     NA                   NA
#> s162                          NA                     NA                   NA
#> s167                          NA                     NA                   NA
#> s173                          NA            -0.02529155                   NA
#> s174                          NA                     NA                   NA
#> s176                          NA            -2.98461811                   NA
#> s177                          NA            -2.98461811                   NA
#> s179                          NA                     NA                   NA
#> s180                          NA                     NA                   NA
#> s181                          NA                     NA                   NA
#> s188                          NA                     NA           0.04104979
#> s197                          NA                     NA                   NA
#> s198                          NA                     NA                   NA
#> s199                          NA                     NA                   NA
#> s212                          NA                     NA                   NA
#> s219                          NA                     NA                   NA
#> s220                          NA                     NA                   NA
#> s221                          NA                     NA                   NA
#> s234                          NA                     NA                   NA
#> s237                          NA                     NA                   NA
#> s240                          NA                     NA                   NA
#> s244                          NA                     NA                   NA
#> s246                 -0.35240174                     NA          -3.03112463
#> s265                          NA                     NA                   NA
#> s270                 -6.66483810            -6.43261730          -5.97381227
#> s271                 -6.66483810            -6.43261730          -5.97381227
#> s274                          NA                     NA                   NA
#> s278                          NA                     NA                   NA
#> s283                          NA                     NA                   NA
#> s287                          NA                     NA                   NA
#> s299                 -3.68356880                     NA          -3.98612399
#> s302                          NA                     NA                   NA
#> s304                          NA                     NA                   NA
#> s307                          NA                     NA                   NA
#> s311                          NA                     NA                   NA
#> s313                          NA                     NA          -0.99066429
#> s322                          NA                     NA                   NA
#> s325                          NA                     NA                   NA
#> s330                          NA                     NA                   NA
#> s341                          NA                     NA                   NA
#> s347                          NA                     NA                   NA
#> s348                          NA                     NA                   NA
#> s356                          NA                     NA                   NA
#> s362                          NA                     NA           0.95245543
#> s374                          NA                     NA                   NA
#> s376                          NA                     NA                   NA
#> s377                          NA                     NA                   NA
#> s381                  0.03161015            -1.03772016                   NA
#> s382                          NA                     NA                   NA
#> s393                 -1.58066105            -2.35553973          -2.23757754
#> s395                          NA                     NA                   NA
#> s396                          NA                     NA                   NA
#> s401                          NA                     NA                   NA
#> s416                          NA            -6.05459464                   NA
#> s418                          NA                     NA                   NA
#> s430                          NA                     NA                   NA
#> s436                 -5.94319031                     NA                   NA
#> s444                          NA                     NA                   NA
#> s448                 -6.98898512            -6.96448577          -7.83055829
#> s451                          NA                     NA                   NA
#> s452                          NA                     NA                   NA
#> s454                          NA                     NA                   NA
#> s457                          NA                     NA                   NA
#> s470                          NA                     NA                   NA
#> s472                          NA                     NA                   NA
#> s474                  1.92784628            -0.90649324           1.19740090
#> s482                          NA                     NA           8.69336266
#> s488                          NA                     NA          -0.84858875
#> s499                          NA                     NA                   NA
#> s500                          NA                     NA                   NA
#>      Phospho_EGF_100min_rep1 Phospho_HGF_10min_rep1 Phospho_HGF_24h_rep1
#> s13                       NA              -4.089779                   NA
#> s16                       NA                     NA                   NA
#> s36                       NA                     NA                   NA
#> s37                       NA                     NA                   NA
#> s40                       NA                     NA                   NA
#> s57                       NA              -1.603229                   NA
#> s58                       NA              -1.603229                   NA
#> s59                       NA              -1.603229                   NA
#> s60                       NA              -1.603229                   NA
#> s64                       NA                     NA                   NA
#> s65                       NA                     NA                   NA
#> s79                       NA                     NA                   NA
#> s83                       NA                     NA                   NA
#> s93                       NA                     NA                   NA
#> s94                       NA                     NA                   NA
#> s114                      NA                     NA                   NA
#> s128                      NA                     NA                   NA
#> s129                      NA                     NA                   NA
#> s136                      NA                     NA                   NA
#> s142                      NA                     NA                   NA
#> s144                      NA                     NA                   NA
#> s145                      NA                     NA                   NA
#> s147                      NA                     NA                   NA
#> s150                      NA                     NA                   NA
#> s151                      NA                     NA                   NA
#> s153                      NA                     NA                   NA
#> s154                      NA                     NA                   NA
#> s155                      NA                     NA                   NA
#> s158                      NA                     NA                   NA
#> s162                      NA                     NA                   NA
#> s167                      NA                     NA                   NA
#> s173                      NA                     NA                   NA
#> s174                      NA                     NA                   NA
#> s176                      NA                     NA                   NA
#> s177                      NA                     NA                   NA
#> s179                      NA                     NA                   NA
#> s180                      NA                     NA                   NA
#> s181                      NA                     NA                   NA
#> s188                      NA                     NA                   NA
#> s197                      NA              -5.310601                   NA
#> s198                      NA                     NA                   NA
#> s199                      NA                     NA                   NA
#> s212                      NA                     NA                   NA
#> s219                      NA                     NA                   NA
#> s220                      NA                     NA                   NA
#> s221                      NA                     NA                   NA
#> s234                      NA                     NA                   NA
#> s237                      NA                     NA                   NA
#> s240                      NA                     NA                   NA
#> s244                      NA                     NA                   NA
#> s246                      NA                     NA                   NA
#> s265                      NA                     NA                   NA
#> s270              -7.9953481              -5.295009                   NA
#> s271              -7.9953481              -5.295009                   NA
#> s274                      NA                     NA           -6.2659013
#> s278                      NA                     NA                   NA
#> s283                      NA                     NA                   NA
#> s287                      NA                     NA                   NA
#> s299              -4.0211310              -3.260359           -2.5015333
#> s302                      NA                     NA                   NA
#> s304                      NA                     NA                   NA
#> s307                      NA                     NA                   NA
#> s311                      NA                     NA                   NA
#> s313                      NA              -2.469246                   NA
#> s322                      NA                     NA                   NA
#> s325                      NA                     NA                   NA
#> s330                      NA                     NA                   NA
#> s341                      NA                     NA                   NA
#> s347                      NA                     NA                   NA
#> s348                      NA                     NA                   NA
#> s356                      NA                     NA                   NA
#> s362                      NA                     NA                   NA
#> s374                      NA                     NA                   NA
#> s376                      NA                     NA                   NA
#> s377                      NA                     NA                   NA
#> s381                      NA              -1.616862           -1.1535401
#> s382                      NA                     NA                   NA
#> s393              -1.0621359              -2.535146           -1.3683106
#> s395                      NA                     NA                   NA
#> s396                      NA                     NA                   NA
#> s401                      NA                     NA                   NA
#> s416              -9.3688472                     NA                   NA
#> s418                      NA                     NA                   NA
#> s430                      NA                     NA                   NA
#> s436                      NA              -7.051126                   NA
#> s444                      NA                     NA                   NA
#> s448              -8.5337133              -8.162944                   NA
#> s451                      NA                     NA                   NA
#> s452               2.9492475                     NA                   NA
#> s454                      NA                     NA                   NA
#> s457                      NA                     NA                   NA
#> s470                      NA                     NA           -0.8770891
#> s472                      NA                     NA                   NA
#> s474               0.3295746               0.388710            0.5693116
#> s482                      NA                     NA                   NA
#> s488                      NA                     NA                   NA
#> s499                      NA                     NA                   NA
#> s500                      NA                     NA                   NA
#>      Phospho_HGF_100min_rep1
#> s13                -3.744997
#> s16                       NA
#> s36                       NA
#> s37                       NA
#> s40                       NA
#> s57                       NA
#> s58                       NA
#> s59                       NA
#> s60                       NA
#> s64                       NA
#> s65                       NA
#> s79                       NA
#> s83                       NA
#> s93                       NA
#> s94                       NA
#> s114                      NA
#> s128                      NA
#> s129                      NA
#> s136                      NA
#> s142                      NA
#> s144                      NA
#> s145                      NA
#> s147                      NA
#> s150                      NA
#> s151                      NA
#> s153                      NA
#> s154                      NA
#> s155                      NA
#> s158                      NA
#> s162                      NA
#> s167                      NA
#> s173                      NA
#> s174                      NA
#> s176                      NA
#> s177                      NA
#> s179                      NA
#> s180                      NA
#> s181                      NA
#> s188                      NA
#> s197                      NA
#> s198                      NA
#> s199                      NA
#> s212                      NA
#> s219                      NA
#> s220                      NA
#> s221                      NA
#> s234                      NA
#> s237                      NA
#> s240                      NA
#> s244                      NA
#> s246                      NA
#> s265                      NA
#> s270              -11.388234
#> s271              -11.388234
#> s274                      NA
#> s278                      NA
#> s283                      NA
#> s287                      NA
#> s299               -1.938292
#> s302                      NA
#> s304                      NA
#> s307                      NA
#> s311                      NA
#> s313                      NA
#> s322                      NA
#> s325                      NA
#> s330                      NA
#> s341                      NA
#> s347                      NA
#> s348                      NA
#> s356                      NA
#> s362                      NA
#> s374                      NA
#> s376                      NA
#> s377                      NA
#> s381                      NA
#> s382                      NA
#> s393               -1.916995
#> s395                      NA
#> s396                      NA
#> s401                      NA
#> s416                      NA
#> s418                      NA
#> s430                      NA
#> s436               -7.965392
#> s444                      NA
#> s448               -9.416012
#> s451                      NA
#> s452                      NA
#> s454                      NA
#> s457                      NA
#> s470                      NA
#> s472                      NA
#> s474               -2.251489
#> s482                      NA
#> s488                      NA
#> s499                      NA
#> s500                      NA
```
