ip_logLike_sub0 = zeros(1000,1);
ip_logLike_sub0 = [-12.7777
-11.6046
-11.7501
-10.3571
-8.23629
-8.09122
-7.8758
-8.11014
-9.56532
-9.35846
-11.2069
-10.9891
-9.61474
-9.32974
-9.08078
-7.43322
-7.54675
-6.70453
-6.98426
-5.12856
-3.84561
-3.38352
-2.44737
-2.11098
-1.77449
-1.81429
-1.87857
-1.79398
-1.57769
-1.63763
-2.12339
-1.74025
-2.00468
-2.08082
-2.19479
-2.15566
-2.23858
-2.81245
-2.53826
-2.35197
-1.9881
-1.57114
-1.42441
-1.31111
-1.00708
-0.989538
-0.967046
-1.27168
-1.52094
-1.53732
-1.18777
-0.610765
-0.649146
-0.522296
-0.583818
-0.425296
-0.402258
-0.454829
-0.43955
-0.548223
-0.543083
-0.508184
-0.452122
-0.408407
-0.393593
-0.360827
-0.33513
-0.354248
-0.227686
-0.216369
-0.249398
-0.313541
-0.260747
-0.187333
-0.217156
-0.188994
-0.123889
-0.119128
-0.117997
-0.112896
-0.149622
-0.132329
-0.168749
-0.226001
-0.298152
-0.363674
-0.365823
-0.309718
-0.380725
-0.361243
-0.290918
-0.292557
-0.228976
-0.304341
-0.219778
-0.15479
-0.167762
-0.100806
-0.126286
-0.0997241
-0.0703151
-0.0572475
-0.0539408
-0.0877266
-0.108571
-0.0887158
-0.0971281
-0.0957689
-0.118177
-0.130946
-0.183204
-0.100824
-0.21158
-0.28153
-0.22362
-0.222057
-0.152774
-0.141897
-0.108776
-0.0915334
-0.107544
-0.0636042
-0.0430591
-0.0511628
-0.0507071
-0.0963464
-0.0770561
-0.0522063
-0.0842501
-0.0845191
-0.0911445
-0.0541555
-0.0594269
-0.0451301
-0.0331649
-0.041016
-0.0666364
-0.0357184
-0.0405311
-0.0303984
-0.0211796
-0.0177139
-0.0119692
-0.013911
-0.0133952
-0.0108508
-0.0155316
-0.00986152
-0.0180774
-0.0303997
-0.0176252
-0.0228285
-0.0154981
-0.0199624
-0.0140361
-0.0152816
-0.0167576
-0.0150416
-0.0118653
-0.0114852
-0.0119963
-0.0112081
-0.0117873
-0.0064504
-0.00672383
-0.0056057
-0.00371279
-0.00267436
-0.00288013
-0.00419048
-0.002713
-0.00176833
-0.00176963
-0.00185719
-0.00218986
-0.00234112
-0.00238666
-0.00219104
-0.00263979
-0.00334675
-0.00297659
-0.00387007
-0.0039116
-0.00365476
-0.00308526
-0.00467469
-0.00334204
-0.00373949
-0.00471731
-0.00464991
-0.00277044
-0.00242412
-0.00237764
-0.00268655
-0.00211668
-0.00191924
-0.0023913
-0.00208094
-0.00178493
-0.00231203
-0.00218373
-0.00296858
-0.00277942
-0.00436471
-0.00211098
-0.00164858
-0.00264592
-0.0018396
-0.00233318
-0.00190092
-0.00270587
-0.00505547
-0.002315
-0.00211819
-0.0023545
-0.00276398
-0.00129388
-0.000812661
-0.000837344
-0.000821246
-0.00116676
-0.000930998
-0.000850943
-0.000718802
-0.000964254
-0.000975666
-0.000824724
-0.000878332
-0.00102289
-0.000949082
-0.000669151
-0.000478988
-0.000565797
-0.000471083
-0.000423141
-0.000581638
-0.000460362
-0.000536482
-0.000601351
-0.000839302
-0.00059007
-0.000344738
-0.000781492
-0.000710537
-0.000515461
-0.000467819
-0.000266228
-0.000354197
-0.000369528
-0.000351798
-0.000371098
-0.000647431
-0.000945728
-0.000877097
-0.00127149
-0.00101019
-0.00139615
-0.00185694
-0.00118379
-0.00123088
-0.00114769
-0.00106743
-0.000884042
-0.00074437
-0.000856577
-0.000723393
-0.000786421
-0.000867549
-0.000780709
-0.000828409
-0.0012632
-0.000869522
-0.000956634
-0.00144855
-0.00097463
-0.00119823
-0.000540461
-0.000518366
-0.000441061
-0.000509646
-0.000544436
-0.000439666
-0.000388659
-0.000481923
-0.000558425
-0.000501757
-0.000413598
-0.000382694
-0.000278614
-0.000263561
-0.000283293
-0.00039763
-0.000583219
-0.000406618
-0.000258444
-0.0003575
-0.00049863
-0.000791292
-0.000758351
-0.000654546
-0.000829026
-0.00104226
-0.000634245
-0.00114235
-0.00098144
-0.00100387
-0.00046657
-0.000334385
-0.000250444
-0.000367935
-0.000250964
-0.000431453
-0.00036172
-0.000408528
-0.000275168
-0.000420251
-0.000477864
-0.000423405
-0.000461688
-0.000717979
-0.00154161
-0.00170098
-0.00169478
-0.00191346
-0.00149065
-0.00133125
-0.00134241
-0.00163385
-0.00186731
-0.0013968
-0.00142891
-0.00114384
-0.00115909
-0.00129846
-0.00167325
-0.00211281
-0.00231895
-0.00295742
-0.00304271
-0.00431773
-0.005861
-0.00531542
-0.00774369
-0.00995489
-0.0115776
-0.0118756
-0.00766713
-0.00750749
-0.0068254
-0.0112502
-0.00802977
-0.0119631
-0.0169219
-0.0204632
-0.02299
-0.0144208
-0.00908582
-0.00635085
-0.00568943
-0.00835042
-0.00896785
-0.0120589
-0.0149791
-0.022836
-0.0197819
-0.0133726
-0.0141447
-0.0162564
-0.0129293
-0.0120269
-0.0123662
-0.00692099
-0.00937999
-0.0131783
-0.0128232
-0.00809419
-0.00825719
-0.0086778
-0.00859291
-0.0112338
-0.0237184
-0.027096
-0.0385678
-0.0252847
-0.01402
-0.0107186
-0.0243016
-0.0173149
-0.0171801
-0.0113436
-0.0155415
-0.0190974
-0.0287066
-0.0481504
-0.0383906
-0.0412277
-0.0409247
-0.0429234
-0.0387596
-0.036372
-0.0368244
-0.0273969
-0.0194702
-0.0221116
-0.0233752
-0.0292806
-0.0329701
-0.0365424
-0.0621992
-0.0688644
-0.170706
-0.152071
-0.20565
-0.265552
-0.239193
-0.290878
-0.249037
-0.207964
-0.210963
-0.226961
-0.168639
-0.248649
-0.328717
-0.365353
-0.253161
-0.281102
-0.240506
-0.264408
-0.126158
-0.116011
-0.108884
-0.0937764
-0.117037
-0.0697599
-0.0875541
-0.0980195
-0.125466
-0.18789
-0.221934
-0.179477
-0.208798
-0.169405
-0.232655
-0.212946
-0.242506
-0.277088
-0.299454
-0.260318
-0.263453
-0.320046
-0.357034
-0.466072
-0.466935
-0.335184
-0.411062
-0.475622
-0.363585
-0.291533
-0.331288
-0.318495
-0.246476
-0.238723
-0.217881
-0.329584
-0.260543
-0.251328
-0.246943
-0.257505
-0.285847
-0.316233
-0.234348
-0.222029
-0.127191
-0.117881
-0.129524
-0.0991105
-0.117376
-0.0973264
-0.139048
-0.1601
-0.187031
-0.246164
-0.316541
-0.261962
-0.300324
-0.387964
-0.329499
-0.323884
-0.364516
-0.362958
-0.278458
-0.343664
-0.443658
-0.427366
-0.480472
-0.436986
-0.696307
-0.76145
-0.897422
-1.06161
-0.998888
-0.919334
-0.937815
-1.05526
-0.803961
-0.887984
-1.02735
-1.24818
-1.08845
-1.33706
-0.968519
-0.967811
-1.1042
-1.10003
-1.18705
-1.0732
-0.98855
-0.787997
-0.994092
-0.862851
-0.98337
-0.735505
-0.553066
-0.438849
-0.526569
-0.630319
-0.614128
-0.580542
-0.768675
-0.455017
-0.444025
-0.567903
-0.609155
-0.809774
-0.646853
-0.571709
-0.443938
-0.353144
-0.332691
-0.345298
-0.505859
-0.562386
-0.566926
-0.576455
-0.566072
-0.656455
-0.927081
-1.06159
-0.816911
-0.582393
-0.68004
-0.708194
-0.682606
-0.625791
-0.439602
-0.48947
-0.33428
-0.326293
-0.410291
-0.52185
-0.558402
-0.811723
-0.739478
-0.678634
-0.799184
-0.746935
-0.698497
-0.572221
-0.746347
-0.56575
-0.597358
-0.497513
-0.486832
-0.386358
-0.414357
-0.326444
-0.464501
-0.54106
-0.342044
-0.389093
-0.416038
-0.262912
-0.305668
-0.303026
-0.264425
-0.166919
-0.226801
-0.243725
-0.203641
-0.201737
-0.27605
-0.385879
-0.341839
-0.281923
-0.354267
-0.277441
-0.200092
-0.291574
-0.25702
-0.180394
-0.199196
-0.206989
-0.231993
-0.210917
-0.165519
-0.173683
-0.290783
-0.436499
-0.373081
-0.296106
-0.323835
-0.346208
-0.314187
-0.30765
-0.23403
-0.251875
-0.19485
-0.149052
-0.0978169
-0.0882122
-0.0925642
-0.0624432
-0.0515474
-0.0622852
-0.0735446
-0.0866996
-0.072758
-0.0465003
-0.0315439
-0.0516668
-0.0626565
-0.0601711
-0.0479203
-0.0816674
-0.138622
-0.145144
-0.121733
-0.0900297
-0.0878285
-0.0995834
-0.131466
-0.140668
-0.11086
-0.0785179
-0.068101
-0.0697007
-0.0655638
-0.0693347
-0.0681581
-0.0529583
-0.051286
-0.0437275
-0.069199
-0.0859782
-0.146767
-0.106976
-0.0766721
-0.0703208
-0.0576499
-0.0467831
-0.0436049
-0.0391786
-0.0346268
-0.0465479
-0.0403837
-0.0537033
-0.0519436
-0.0443206
-0.0244573
-0.0492042
-0.0658671
-0.0735037
-0.0705807
-0.0939645
-0.0835676
-0.107023
-0.104238
-0.144011
-0.0902601
-0.0783478
-0.0406134
-0.0325281
-0.0409762
-0.0403881
-0.0490376
-0.0488355
-0.0591429
-0.0641128
-0.0973548
-0.116751
-0.0781301
-0.107514
-0.14589
-0.117813
-0.0932077
-0.114331
-0.14729
-0.208983
-0.13802
-0.225543
-0.26564
-0.390738
-0.424839
-0.293905
-0.205622
-0.2121
-0.186237
-0.259394
-0.317133
-0.30485
-0.346734
-0.333385
-0.217142
-0.194595
-0.138536
-0.151165
-0.139758
-0.167261
-0.165487
-0.176066
-0.137447
-0.153938
-0.173118
-0.202173
-0.191983
-0.194258
-0.189334
-0.14595
-0.144248
-0.232651
-0.327983
-0.337804
-0.412403
-0.529835
-0.327036
-0.200693
-0.297553
-0.350465
-0.363497
-0.321306
-0.220932
-0.333938
-0.26903
-0.348989
-0.34393
-0.388136
-0.615551
-0.562522
-0.631763
-0.643195
-0.525913
-0.499283
-0.468849
-0.402532
-0.368955
-0.24573
-0.351406
-0.376715
-0.426938
-0.432963
-0.43429
-0.445924
-0.472039
-0.363669
-0.48899
-0.356722
-0.401529
-0.561686
-0.475113
-0.436465
-0.528915
-0.767686
-0.605124
-0.449498
-0.449478
-0.444717
-0.534106
-0.601225
-0.529767
-0.407527
-0.453401
-0.414233
-0.472311
-0.241419
-0.203031
-0.18324
-0.164557
-0.124454
-0.169545
-0.219278
-0.214624
-0.210663
-0.143972
-0.112039
-0.0824837
-0.0721919
-0.0544802
-0.0314305
-0.049733
-0.0391909
-0.0727568
-0.0758572
-0.0546142
-0.0326329
-0.0508677
-0.0362638
-0.0344425
-0.0441607
-0.0311406
-0.0388539
-0.0304954
-0.0291654
-0.0444285
-0.0594642
-0.0618352
-0.0631732
-0.04037
-0.0357515
-0.0410133
-0.0492045
-0.0795971
-0.0892516
-0.0945775
-0.0769741
-0.0684092
-0.0749242
-0.0813966
-0.107982
-0.112301
-0.0955155
-0.125571
-0.0700127
-0.058579
-0.0636483
-0.0604464
-0.0422733
-0.0461965
-0.028748
-0.0322562
-0.0282026
-0.0482324
-0.0410643
-0.0385618
-0.036695
-0.0347015
-0.0388533
-0.0367615
-0.0315213
-0.0326974
-0.0231407
-0.034384
-0.030542
-0.0324475
-0.0297138
-0.0371185
-0.0354464
-0.0392077
-0.0575489
-0.0353347
-0.0492125
-0.063564
-0.0544458
-0.0321958
-0.0318705
-0.0328148
-0.0274903
-0.0439678
-0.0329328
-0.0289937
-0.0129382
-0.0129604
-0.0215552
-0.01786
-0.0133732
-0.0101464
-0.0156432
-0.0176348
-0.0306626
-0.0391267
-0.0206692
-0.0254832
-0.0190403
-0.0182596
-0.0125682
-0.0120605
-0.0133865
-0.0101477
-0.00890221
-0.00646918
-0.00905097
-0.0104232
-0.00984818
-0.0104524
-0.0126049
-0.018433
-0.0306661
-0.0472015
-0.0525175
-0.0878078
-0.0712712
-0.0980334
-0.0783553
-0.0888987
-0.13842
-0.136981
-0.219297
-0.240958
-0.279542
-0.179977
-0.23666
-0.230762
-0.226584
-0.155104
-0.137668
-0.180526
-0.205714
-0.128731
-0.133427
-0.151047
-0.140167
-0.156641
-0.148202
-0.202717
-0.208723
-0.205682
-0.261048
-0.160079
-0.177979
-0.0650631
-0.0590674
-0.0345846
-0.0394799
-0.0316373
-0.0502322
-0.0418496
-0.0516104
-0.0367741
-0.0376744
-0.0462431
-0.0403613
-0.0558155
-0.0641444
-0.0596425
-0.0381705
-0.0503696
-0.0357267
-0.0241232
-0.0298023
-0.0317938
-0.0252184
-0.0263757
-0.0239122
-0.0155961
-0.0116384
-0.0151513
-0.0165704
-0.019064
-0.0220625
-0.018951
-0.0195189
-0.014367
-0.0197581
-0.0142828
-0.0130118
-0.0141181
-0.0119893
-0.0142476
-0.00980608
-0.0143111
-0.0148495
-0.00999828
-0.0108101
-0.00880264
-0.00832721
-0.0130228
-0.0144827
-0.0121936
-0.0120165
-0.01537
-0.0103502
-0.0096474
-0.00963634
-0.00681758
-0.00728924
-0.00772564
-0.0155739
-0.0130365
-0.0102608
-0.0097531
-0.0152266
-0.0219038
-0.0211501
-0.023846
-0.017059
-0.0204494
-0.0404708
-0.0468004
-0.0512725
-0.0707579
-0.145436
];
ip_logLike_sub0 = zeros(1000,1);
ip_logLike_sub0 = [4.99131
4.89994
4.82361
4.98226
5.23732
5.41543
5.32171
5.29616
5.24881
4.67942
4.95567
4.63797
4.88563
4.52402
4.76326
4.99881
5.04953
5.03687
5.03976
5.33408
5.09469
4.78193
4.78507
4.5936
4.4191
4.06195
4.01389
4.4262
4.46105
5.02763
5.07783
5.48723
5.7834
5.72755
5.47593
5.37449
5.64066
5.76436
5.57192
6.02328
5.94014
5.85248
5.68246
5.74175
5.69139
5.43741
5.09722
5.07979
4.90763
5.05992
5.43532
5.24875
4.89415
4.93242
5.27251
5.59111
5.81886
5.71494
5.63727
5.63924
6.09114
6.26633
6.27793
5.96418
5.92085
5.58464
5.43512
5.44756
5.19684
5.22751
5.07946
4.9303
5.01968
5.22541
5.43709
5.18034
5.11827
5.19281
5.33871
5.53827
4.98368
5.07031
4.89459
4.92021
4.93142
5.09562
5.16472
5.19827
5.4815
5.57981
5.54625
5.64769
5.61965
5.50413
5.54544
5.81993
5.63005
5.49293
5.761
5.66741
5.47451
5.14344
5.40375
5.5085
5.58723
5.29602
5.50182
5.20774
5.19694
5.24422
5.52449
5.56017
5.01995
5.05143
5.04143
5.06987
4.77583
4.68137
4.52018
4.34484
4.15178
4.40365
4.64379
4.57225
4.63696
4.75877
4.75777
5.16428
5.4463
5.39167
5.25248
5.32676
4.95675
4.80261
4.49003
4.40556
4.54706
4.41417
4.45552
4.25634
4.25447
4.03709
4.06362
3.93501
3.82841
3.94091
4.11932
4.26986
4.42448
4.4904
4.58269
4.69826
4.47379
4.31404
4.27809
3.91902
3.8555
3.8717
3.66775
3.57709
3.72219
3.51998
3.8067
3.79766
3.45706
3.46862
3.62403
3.56137
3.43544
3.69151
3.8707
3.47559
3.23291
3.32456
3.38336
3.56756
3.46784
3.55032
3.46209
3.62051
3.64845
4.00084
4.08024
4.09523
3.77193
3.79306
4.20724
4.11437
4.54972
4.43367
4.69933
4.40038
4.29029
4.42455
4.26611
4.04047
4.44666
4.20629
4.17243
4.02905
4.40063
4.80309
4.51086
4.95547
4.58643
4.33139
4.46647
4.38684
4.56861
4.68213
4.57911
5.03365
4.70965
4.73757
4.6592
4.86007
4.36007
4.03012
3.82054
3.93317
4.11604
3.97923
4.0024
3.87749
4.0567
4.05103
3.93371
3.98452
3.97975
4.06458
3.86255
3.56932
3.72298
3.58708
3.5074
3.74344
3.56999
3.68353
3.76813
4.01439
3.7541
3.35527
3.96181
3.89155
3.65389
3.42299
3.06228
3.2716
3.33755
3.18258
3.24285
3.73231
4.0931
4.148
4.09636
4.13005
4.33613
4.12563
4.30302
4.2195
4.35766
3.98481
3.85352
3.85106
3.88486
3.91275
4.24109
4.18242
4.02093
4.08722
4.3179
4.21689
4.53958
4.66239
4.49217
4.2156
4.48655
4.40466
4.28535
4.44083
4.31704
4.51763
4.41794
4.59363
4.72455
4.41606
4.37805
3.9897
4.11086
4.05431
4.39726
4.72699
4.46317
4.26077
4.20159
4.66729
4.86494
5.07591
4.93331
5.03754
5.10358
4.98843
5.20027
4.95345
5.29471
4.71086
4.26077
3.97388
4.12845
4.17538
4.21763
4.19746
4.44003
4.21246
4.47592
4.57128
4.58561
4.76761
5.19512
5.32435
5.56792
5.55596
5.62867
5.68287
5.60204
5.58281
5.74587
5.72554
5.59426
5.61156
5.49054
5.51195
5.4176
5.58497
5.7654
5.96233
5.97593
6.05998
6.38864
6.47588
6.35453
6.06037
6.3548
6.40775
6.53595
6.48299
6.62787
6.45866
6.5797
6.49981
6.40757
6.69318
6.75368
6.59869
6.72934
6.98776
6.83973
6.81694
6.95096
7.16255
7.07401
6.88901
6.934
6.83356
7.12113
7.20257
7.14987
7.04454
6.86379
6.73467
6.84941
6.95387
7.03106
7.00265
7.1015
7.16652
6.88476
7.14573
6.89242
6.96824
6.90353
6.79396
6.97696
6.90596
7.12766
7.38974
7.51928
7.62776
7.69656
7.87069
8.02758
8.03842
7.94773
8.1391
8.20901
8.16094
7.98773
8.138
8.27343
8.32958
8.3718
8.41937
8.50569
8.60557
8.57058
8.48677
8.44982
8.41908
8.2744
8.41183
8.35215
8.39006
8.55575
8.63816
8.76678
8.7989
8.72687
8.65118
8.6562
8.62367
8.64075
8.63983
8.56987
8.64326
8.56603
8.61435
8.47868
8.49417
8.35948
8.34865
8.45693
8.50676
8.47892
8.43325
8.67994
8.76053
8.66381
8.67318
8.66832
8.66095
8.60475
8.64514
8.65152
8.47909
8.50566
8.484
8.54507
8.6155
8.65669
8.73592
8.73083
8.76068
8.85955
8.83878
8.80699
8.89353
8.72001
8.70538
8.66546
8.73512
8.75026
8.64751
8.56256
8.4832
8.47802
8.47489
8.46228
8.37031
8.42572
8.38666
8.21851
8.04017
7.93149
8.02641
8.05708
8.09212
7.94999
8.13015
7.8797
7.66147
7.38468
7.34886
7.40734
7.2131
7.40326
7.44622
7.4482
7.45078
7.25508
7.08328
7.11911
7.24313
7.25889
7.08078
6.97009
7.0829
7.142
7.34954
7.44542
7.17882
6.87249
6.85891
6.87826
6.68712
6.85472
7.06313
7.09737
7.14905
7.10875
7.30225
7.22481
7.3734
7.33204
7.17743
7.11927
6.67017
6.66302
6.70785
7.07737
7.03161
7.11572
7.02642
6.9918
7.33157
7.40773
7.29941
7.34893
7.14327
7.30724
7.14952
7.28764
7.25578
7.4003
7.48053
7.61009
7.6539
7.53376
7.37404
7.5981
7.65523
7.69282
7.75713
7.89268
7.74882
7.72873
7.74878
7.87113
7.6045
7.69039
7.76827
7.704
7.79545
7.72699
7.77782
7.83591
7.82231
7.53019
7.45869
7.58952
7.75318
7.82358
7.73074
7.68949
7.71459
7.59739
7.67726
7.64263
7.45574
7.54789
7.3847
7.08529
7.01209
7.1697
7.11401
7.00523
6.94123
7.01146
7.16387
7.39276
7.05956
6.89006
6.77734
6.77721
6.83007
6.83102
6.99424
7.14897
6.91808
6.62654
6.51008
6.70779
6.76784
6.73851
6.53433
6.34085
6.28715
6.12463
5.95199
6.24045
6.30867
5.82415
5.86356
5.83094
5.44788
5.39736
5.55342
5.76388
6.13845
6.28584
6.17202
6.5653
6.4494
6.42884
6.49108
6.42704
6.21131
6.17129
5.88198
6.06684
6.17775
5.82386
5.6053
6.11491
6.19143
6.49964
6.66739
6.89579
7.04312
7.17114
7.42448
7.34113
7.32092
7.39016
7.5377
7.77859
7.78898
7.84868
7.77375
7.78488
7.92253
7.84073
7.89427
7.63034
7.61663
7.58919
7.56807
7.6815
7.76226
7.95434
7.87379
7.66415
7.54653
7.51642
7.52727
7.64996
7.5919
7.71345
7.87256
7.74459
7.33463
7.3775
7.59487
7.61595
7.44239
7.27375
7.31758
7.31451
7.57019
7.29273
7.37704
7.30905
7.34906
7.51561
7.41856
7.57777
7.50105
7.59657
7.70043
7.80081
7.78888
7.89864
7.78222
7.79188
7.75358
7.58646
7.51475
7.61194
7.70579
7.81831
7.83876
7.66417
7.74981
7.95304
8.09558
8.1939
8.13464
7.95154
7.87954
8.0017
7.91177
8.05777
8.11098
8.2046
8.22158
8.42309
8.49436
8.53907
8.46528
8.59097
8.65559
8.53937
8.53643
8.5509
8.54197
8.34299
8.42191
8.46511
8.37162
8.42166
8.41794
8.45512
8.46962
8.50765
8.60589
8.50102
8.48528
8.44211
8.39493
8.48524
8.42257
8.50127
8.51668
8.6092
8.54335
8.49858
8.44645
8.521
8.66235
8.72005
8.81444
8.90448
8.94759
8.93019
8.84887
8.93796
8.83942
8.98233
8.9451
8.85703
8.95114
8.9593
8.97959
8.93617
8.96293
8.83786
8.76997
8.88657
8.76152
8.79649
8.85368
8.84409
8.8484
8.97427
8.92555
8.85971
8.78627
8.77181
8.80072
8.65631
8.71329
8.75812
8.78924
8.82542
8.87431
8.91614
8.98413
8.92819
8.9516
8.94302
8.88387
8.81908
8.83317
8.69365
8.62775
8.42621
8.3693
8.42504
8.42204
8.29993
8.40359
8.61592
8.6195
8.59351
8.78984
8.60202
8.4645
8.53808
8.38855
8.3559
8.4135
8.52085
8.59663
8.74228
8.9049
8.72014
8.77134
8.81957
8.73135
8.76351
8.8489
8.89167
9.01082
9.04247
9.1089
9.14408
9.13851
9.19468
9.14543
9.14698
9.14601
9.16757
9.30203
9.27098
9.29057
9.25561
9.24335
9.24831
9.30337
9.33207
9.30544
9.30057
9.35568
9.25403
9.19548
9.26532
9.23429
9.14496
9.18608
9.11328
9.11761
9.10716
9.2033
9.17958
9.12084
9.14515
9.11177
9.12154
9.13889
9.09184
9.11323
9.09304
9.16738
9.14563
9.15679
9.14052
9.18121
9.1729
9.191
9.25737
9.17233
9.23078
9.27393
9.24802
9.15536
9.15349
9.15885
9.12596
9.21123
9.15951
9.13595
8.97483
8.98873
9.09147
9.0545
8.99532
8.94249
8.90188
8.78659
8.857
8.8301
8.82086
8.77782
8.63879
8.54755
8.37776
8.17063
8.1432
8.10314
8.00181
7.8494
7.90712
7.92635
7.98045
8.0545
8.10414
8.06372
7.94838
7.83104
7.86942
8.111
8.22395
8.06761
8.06393
7.97223
7.85949
7.7042
7.84898
7.77017
7.86453
7.77631
7.9484
8.03114
8.00419
7.79437
7.65659
7.50992
7.31242
7.15046
7.33862
7.2036
7.27656
7.31634
7.48597
7.54197
7.75478
7.72164
7.56561
7.84224
7.93611
7.93321
7.95243
7.80156
7.71298
7.69114
7.76229
7.69652
7.8959
7.96942
7.95815
7.97772
7.81906
7.95103
8.02256
8.13504
8.1346
8.00726
8.01357
8.16265
8.1304
8.16061
8.20174
8.04633
8.25533
8.12561
8.10517
8.27454
8.35067
8.45411
8.52635
8.58651
8.57024
8.53246
8.45594
8.47116
8.25065
8.31977
8.40614
8.42727
8.46648
8.62356
8.42423
8.5256
8.451
8.69682
8.74584
8.86717
8.79987
8.79153
8.78325
8.75856
8.74687
8.71534
8.6937
8.76043
8.60891
8.70217
8.76931
8.83089
8.74941
8.83303
8.77774
8.87849
8.95628
8.92985
9.03025
9.16355
9.16465
9.20814
9.26183
9.38175
];
