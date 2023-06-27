library(tinytest)
library(BayesRep)

## formatBF
## -----------------------------------------------------------------------------
bf <- c(1/300.01, 1/3, 1/2, 1/1.01, 1, 1.01, 2, 3, 300.01, NaN, NA)
bfexpected <- c("1/300", "1/3", "1/2", "1", "1", "1", "2", "3", "300", NA, NA)
expect_equal(formatBF(bf), bfexpected, info = "formatBF works correctly")

## BFs and BFr
## -----------------------------------------------------------------------------
## results from Table 1 in Pawel and Held (2022), extracted with dput
tablePaper <- structure(list(to = c(1.14417264341251, 1.96382457556354,
                                    0.818339255457479, 0.582737811761101,
                                    0.226778691550525, 0.696711646322321,
                                    0.488611206901194, 0.485326762453132,
                                    0.817287557101302, 0.203339606989006,
                                    0.743741868967827, 0.285393193419145,
                                    0.396018834302785, 0.141926757107812,
                                    0.27657049025531, 0.385698718807712,
                                    0.274109959311116, 0.275860006439804,
                                    0.297627922877741, 0.409601469520218,
                                    1.07988209588739),
                             so = c(0.164398987305357, 0.288675134594813,
                                    0.192450089729875, 0.144337567297406,
                                    0.0755928946018454, 0.164398987305357,
                                    0.185695338177052, 0.218217890235992,
                                    0.288675134594813, 0.0712470499879096,
                                    0.129099444873581, 0.120385853085769,
                                    0.174077655955698, 0.054232614454664,
                                    0.140028008402801, 0.123091490979333,
                                    0.137360563948689, 0.10976425998969,
                                    0.136082763487954, 0.164398987305357,
                                    0.242535625036333),
                             tr = c(1.19476271568933, 1.18451560871901,
                                    0.683204796807765, 0.37746393971232,
                                    0.184400291006785, 0.404898618376599,
                                    0.370800610719356, 0.671862945083057,
                                    0.467849160696654, 0.116566798991139,
                                    0.358414852842928, 0.147477389849293,
                                    0.150864396208723, 0.025792251826568,
                                    0.0627861994145097, 0.0496313025904017,
                                    -0.0150068664441846, -0.0272611982545971,
                                    -0.0354850035691282, -0.0463558676741315,
                                    -0.0985013285629837),
                             sr = c(0.229415733870562, 0.301511344577764,
                                    0.166666666666667, 0.127000127000191,
                                    0.0497518595104995, 0.147441956154897,
                                    0.107832773203438, 0.104257207028537,
                                    0.105999788000636, 0.045786854649563,
                                    0.160128153805087, 0.0645497224367903,
                                    0.0574484989621426, 0.021652326748721,
                                    0.040961596025952, 0.0657951694959769,
                                    0.040291148201269, 0.0375029300308675,
                                    0.0435194139889245, 0.0594438298277764,
                                    0.114707866935281),
                             bfs = c(0.0000530477419198753, 0.0127778208784195,
                                     0.0224319045241547, 0.117135677207107,
                                     0.144938815030417, 0.179077965986501,
                                     0.255508503862225, 0.309170196577008,
                                     0.322063749958956, 0.400925337910988,
                                     0.629934964682855, 0.637603992664161,
                                     0.849104156875212, NaN, NaN, NaN, NaN, NaN,
                                     NaN, NaN, NaN),
                             bfr = c(0.00000161275109902623,
                                     0.00352146020839239, 0.000394726943693376,
                                     0.0323155569085302, 0.00211057610414055,
                                     0.082609423824399, 0.00626432422526435,
                                     0.00000000299745740619506,
                                     0.000325684031619899, 0.122358606947869,
                                     0.606604189711098, 0.259071929155332,
                                     0.248184710610032, 9.58539878281878,
                                     3.21906462367827, 28.9728262321521,
                                     25.4825254562925, 72.2003426853402,
                                     35.6721291232169, 65.1067536780313,
                                     24995.5941001993),
                             study = c("Hauser et al. (2014)",
                                       "Aviezer et al. (2012)",
                                       "Wilson et al. (2014)",
                                       "Derex et al. (2013)",
                                       "Gneezy et al. (2014)",
                                       "Karpicke and Blunt (2011)",
                                       "Morewedge et al. (2010)",
                                       "Kovacs et al. (2010)",
                                       "Duncan et al. (2012)",
                                       "Nishi et al. (2015)",
                                       "Janssen et al. (2010)",
                                       "Balafoutas and Sutter (2012)",
                                       "Pyc and Rawson (2010)",
                                       "Rand et al. (2012)",
                                       "Ackerman et al. (2010)",
                                       "Sparrow et al. (2011)",
                                       "Shah et al. (2012)",
                                       "Kidd and Castano (2013)",
                                       "Gervais and Norenzayan (2012)",
                                       "Lee and Schwarz (2010)",
                                       "Ramirez and Beilock (2011)")),
                        row.names = c(NA, -21L), class = "data.frame")
tablePaper$zo <- tablePaper$to/tablePaper$so
tablePaper$zr <- tablePaper$tr/tablePaper$sr
tablePaper$c <- tablePaper$so^2/tablePaper$sr^2

## BFs tests
bfs1 <- BFs(to = tablePaper$to, so = tablePaper$so, tr = tablePaper$tr,
            sr = tablePaper$sr)
expect_equal(bfs1, tablePaper$bfs,
             info = "BFs (computed with to, so, tr, sr) as in Table 1 from Pawel and Held (2022)")
bfs2 <- BFs(zo = tablePaper$zo, zr = tablePaper$zr, c = tablePaper$c)
expect_equal(bfs2, tablePaper$bfs,
             info = "BFs (computed with zo, zr, c) as in Table 1 from Pawel and Held (2022)")

## BFr tests
bfr1 <- BFr(to = tablePaper$to, so = tablePaper$so, tr = tablePaper$tr,
            sr = tablePaper$sr, ss = 0)
expect_equal(bfr1, tablePaper$bfr,
             info = "BFr (computed with to, so, tr, sr, ss) as in Table 1 from Pawel and Held (2022)")
bfr2 <- BFr(zo = tablePaper$zo, zr = tablePaper$zr, c = tablePaper$c, g = 0)
expect_equal(bfr2, tablePaper$bfr,
             info = "BFr (computed with zo, zr, c, g) as in Table 1 from Pawel and Held (2022)")

## BFrlogOR and BFslogOR
## -----------------------------------------------------------------------------
## results from Table 1 in Pawel and Held (2022), extracted with dput
logORtable <- structure(list(study = c("Balafoutas and Sutter (2012), Science",
                                       "Gneezy et al. (2014), Science",
                                       "Hauser et al. (2014), Nature"),
                             ao = c(21, 65, 20),
                             bo = c(15, 26, 0),
                             co = c(11, 43, 4),
                             do = c(25, 44, 16),
                             ar = c(63, 147, 11),
                             br = c(60, 55, 0),
                             cr = c(44, 113, 2),
                             dr = c(76, 92, 9),
                             bfsint = c(0.633250599371189, 0.133503588926685,
                                        NaN),
                             bfshyg = c(0.633250599371083, NaN, 0.000955394729787569),
                             bfrepint = c(0.258716394614373, 0.00181587960505, NaN),
                             bfrephyg = c(0.258716394614278, NaN,
                                          0.0000595783977455087)),
                        class = "data.frame", row.names = c(NA, -3L))

## ## BFslogOR tests (take too long, don't run by default)
## bfslogor1 <- with(logORtable, BFslogOR(ao = ao, bo = bo, co = co, do = do,
##                                        ar = ar, br = br, cr = cr, dr = dr,
##                                        method = "integration"))
## expect_equal(bfslogor1, logORtable$bfsint,
##              info = "BFslogOR (integration) as in Table 1 from Pawel and Held (2022)")
## bfslogor2 <- with(logORtable, BFslogOR(ao = ao, bo = bo, co = co, do = do,
##                                        ar = ar, br = br, cr = cr, dr = dr,
##                                        method = "hypergeo"))
## expect_equal(bfslogor2, logORtable$bfshyg,
##              info = "BFslogOR (hypergeo) as in Table 1 from Pawel and Held (2022)")

## BFrlogOR tests
bfrlogor1 <- with(logORtable, BFrlogOR(ao = ao, bo = bo, co = co, do = do,
                                       ar = ar, br = br, cr = cr, dr = dr,
                                       method = "integration", ss = 0))
expect_equal(bfrlogor1, logORtable$bfrepint,
             info = "BFrlogOR (integration) as in Table 1 from Pawel and Held (2022)")
bfrlogor2 <- with(logORtable, BFrlogOR(ao = ao, bo = bo, co = co, do = do,
                                       ar = ar, br = br, cr = cr, dr = dr,
                                       method = "hypergeo", ss = 0))
expect_equal(bfrlogor2, logORtable$bfrephyg,
             info = "BFrlogOR (hypergeo) as in Table 1 from Pawel and Held (2022)")



## BFrSMD and BFsSMD
## -----------------------------------------------------------------------------
## results from Table 1 in Pawel and Held (2022), extracted with dput
smdTable <- structure(list(study = c("Ackerman et al. (2010), Science",
                                     "Gervais and Norenzayan (2012), Science",
                                     "Karpicke and Blunt (2011), Science",
                                     "Kidd and Castano (2013), Science",
                                     "Lee and Schwarz (2010), Science",
                                     "Morewedge et al. (2010), Science",
                                     "Nishi et al. (2015), Nature",
                                     "Pyc and Rawson (2010), Science",
                                     "Ramirez and Beilock (2011), Science",
                                     "Rand et al. (2012), Nature",
                                     "Shah et al. (2012), Science",
                                     "Wilson et al. (2014), Science",
                                     "Aviezer et al. (2012), Science",
                                     "Duncan et al. (2012), Science",
                                     "Kovacs et al. (2010), Science",
                                     "Sparrow et al. (2011), Science"),
                           type = c("SMD", "SMD", "SMD", "SMD", "SMD", "SMD",
                                    "SMD", "SMD", "SMD", "SMD", "SMD", "SMD",
                                    "SM", "SM", "SM", "SM"),
                           to = c(2.02, 2.24, 4.65, 2.53, 2.6, 2.78, 2.68, 2.37,
                                  5.53, 2.45, 2.04, 4.83, 13.07, 3.41, 2.42,
                                  3.26),
                           n1o = c(26, 26, 20, 43, 21, 16, 10, 18, 10, 175, 26,
                                   15, 15, 15, 24, 69),
                           n2o = c(28, 31, 20, 43, 19, 16, 10, 18, 10, 168, 30,
                                   15, NA, NA, NA, NA),
                           tr = c(1.5351, -0.82, 2.8825, -0.726, -0.78, 3.54,
                                  2.53, 2.64, -0.7352, 1.01, -0.373, 4.49,
                                  5.342, 4.6276, 7.0152, 0.7579),
                           n1r = c(296, 262, 23, 349, 147, 44, 24, 156, 45,
                                   1058, 298, 20, 14, 92, 95, 234),
                           n2r = c(303, 269, 26, 365, 139, 45, 24, 150, 34,
                                   1078, 321, 19, NA, NA, NA, NA),
                           bfsSMD = c(NaN, NaN, 0.201393824351508, NaN, NaN,
                                      0.251915218801748, 0.461389422414924,
                                      0.847327600666465, NaN, NaN, NaN,
                                      0.0284008547136663, 0.098330118343855,
                                      0.320262992627344, 0.26512121011443, NaN),
                           bfrepSMD = c(3.2425436178894, 36.8293701822416,
                                        0.0865954100114667, 68.5206475865596,
                                        68.7768313720538, 0.00639050298948843,
                                        0.131877512393028, 0.25088787491907,
                                        17376.1099486192, 9.72328315095483,
                                        26.1613308415308, 0.000538516055919849,
                                        0.0243407752343321,
                                        0.000414919994401764,
                                        0.00000000744262267585224,
                                        31.8489025486167)),
                      row.names = c(NA, -16L), class = "data.frame")
smd1 <- subset(smdTable, type == "SMD")
smd2 <- subset(smdTable, type == "SM")

## BFsSMD (put some tolerance because optimization has slightly changed)
bfssmd1 <- with(smd1, BFsSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r,
                             n2r = n2r, type = "two.sample"))
expect_equal(log(bfssmd1), log(smd1$bfsSMD), tolerance = 0.001,
             info = "BFsSMD (two.sample) as in Table 1 from Pawel and Held (2022)")
bfssmd2 <- with(smd2, BFsSMD(to = to, no = n1o, tr = tr, nr = n1r,
                             type = "one.sample"))
expect_equal(log(bfssmd2), log(smd2$bfsSMD), tolerance = 0.001,
             info = "BFsSMD (one.sample) as in Table 1 from Pawel and Held (2022)")

## BFrSMD
bfrsmd1 <- with(smd1, BFrSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r,
                             n2r = n2r, type = "two.sample", ss = 0))
expect_equal(bfrsmd1, smd1$bfrepSMD,
             info = "BFrSMD (two.sample) as in Table 1 from Pawel and Held (2022)")
bfrsmd2 <- with(smd2, BFrSMD(to = to, no = n1o, tr = tr, nr = n1r,
                             type = "one.sample", ss = 0))
expect_equal(bfrsmd2, smd2$bfrepSMD,
             info = "BFrSMD (one.sample) as in Table 1 from Pawel and Held (2022)")
