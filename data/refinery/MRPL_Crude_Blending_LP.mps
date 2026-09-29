NAME          MRPL_Crude_Blending_LP
OBJSENSE
  MAX
ROWS
 N  OBJ
 L  CDU_Capacity
 E  Yield_Naphtha
 E  Yield_JetKero
 E  Yield_Gasoil
 L  Gasoil_Sulfur_Spec
 G  Min_Gasoil_Demand
COLUMNS
    Crude_ArabLight OBJ        -72
    Crude_ArabLight CDU_Capacity 1
    Crude_ArabLight Yield_Naphtha 0.22
    Crude_ArabLight Yield_JetKero 0.29999999999999999
    Crude_ArabLight Yield_Gasoil 0.41999999999999998
    Crude_ArabLight Gasoil_Sulfur_Spec 0.041999999999999989
    Crude_BonnyLight OBJ        -78
    Crude_BonnyLight CDU_Capacity 1
    Crude_BonnyLight Yield_Naphtha 0.28000000000000003
    Crude_BonnyLight Yield_JetKero 0.34000000000000002
    Crude_BonnyLight Yield_Gasoil 0.34999999999999998
    Crude_BonnyLight Gasoil_Sulfur_Spec -0.12249999999999998
    Crude_IranHeavy OBJ        -68
    Crude_IranHeavy CDU_Capacity 1
    Crude_IranHeavy Yield_Naphtha 0.17999999999999999
    Crude_IranHeavy Yield_JetKero 0.25
    Crude_IranHeavy Yield_Gasoil 0.47999999999999998
    Crude_IranHeavy Gasoil_Sulfur_Spec 0.192
    Sell_Naphtha OBJ        95
    Sell_Naphtha Yield_Naphtha -1
    Sell_JetKero OBJ        105
    Sell_JetKero Yield_JetKero -1
    Sell_Gasoil OBJ        102
    Sell_Gasoil Yield_Gasoil -1
    Sell_Gasoil Min_Gasoil_Demand 1
RHS
    RHS1       CDU_Capacity 1000
    RHS1       Min_Gasoil_Demand 300
RANGES
BOUNDS
 UP BND1       Crude_ArabLight 500
 UP BND1       Crude_BonnyLight 500
 UP BND1       Crude_IranHeavy 500
 UP BND1       Sell_Naphtha 1000
 UP BND1       Sell_JetKero 1000
 UP BND1       Sell_Gasoil 1000
ENDATA
