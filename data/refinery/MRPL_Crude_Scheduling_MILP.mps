NAME          MRPL_Crude_Scheduling_MILP
OBJSENSE
  MIN
ROWS
 N  OBJ
 E  OneTank_t1
 G  CDU_Demand_t1
 L  MaxRate_k0_t0
 G  MinRate_k0_t0
 L  MaxRate_k1_t0
 G  MinRate_k1_t0
 L  MaxRate_k2_t0
 G  MinRate_k2_t0
 E  OneTank_t2
 G  CDU_Demand_t2
 L  MaxRate_k0_t1
 G  MinRate_k0_t1
 L  MaxRate_k1_t1
 G  MinRate_k1_t1
 L  MaxRate_k2_t1
 G  MinRate_k2_t1
 E  OneTank_t3
 G  CDU_Demand_t3
 L  MaxRate_k0_t2
 G  MinRate_k0_t2
 L  MaxRate_k1_t2
 G  MinRate_k1_t2
 L  MaxRate_k2_t2
 G  MinRate_k2_t2
 E  OneTank_t4
 G  CDU_Demand_t4
 L  MaxRate_k0_t3
 G  MinRate_k0_t3
 L  MaxRate_k1_t3
 G  MinRate_k1_t3
 L  MaxRate_k2_t3
 G  MinRate_k2_t3
 L  TankInv_k1
 L  TankInv_k2
 L  TankInv_k3
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k1_t1 OBJ        50
    FeedBin_k1_t1 OneTank_t1 1
    FeedBin_k1_t1 MaxRate_k0_t0 -120
    FeedBin_k1_t1 MinRate_k0_t0 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k1_t1 OBJ        10
    Flow_k1_t1 CDU_Demand_t1 1
    Flow_k1_t1 MaxRate_k0_t0 1
    Flow_k1_t1 MinRate_k0_t0 1
    Flow_k1_t1 TankInv_k1 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k1_t2 OBJ        50
    FeedBin_k1_t2 OneTank_t2 1
    FeedBin_k1_t2 MaxRate_k0_t1 -120
    FeedBin_k1_t2 MinRate_k0_t1 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k1_t2 OBJ        10
    Flow_k1_t2 CDU_Demand_t2 1
    Flow_k1_t2 MaxRate_k0_t1 1
    Flow_k1_t2 MinRate_k0_t1 1
    Flow_k1_t2 TankInv_k1 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k1_t3 OBJ        50
    FeedBin_k1_t3 OneTank_t3 1
    FeedBin_k1_t3 MaxRate_k0_t2 -120
    FeedBin_k1_t3 MinRate_k0_t2 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k1_t3 OBJ        10
    Flow_k1_t3 CDU_Demand_t3 1
    Flow_k1_t3 MaxRate_k0_t2 1
    Flow_k1_t3 MinRate_k0_t2 1
    Flow_k1_t3 TankInv_k1 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k1_t4 OBJ        50
    FeedBin_k1_t4 OneTank_t4 1
    FeedBin_k1_t4 MaxRate_k0_t3 -120
    FeedBin_k1_t4 MinRate_k0_t3 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k1_t4 OBJ        10
    Flow_k1_t4 CDU_Demand_t4 1
    Flow_k1_t4 MaxRate_k0_t3 1
    Flow_k1_t4 MinRate_k0_t3 1
    Flow_k1_t4 TankInv_k1 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k2_t1 OBJ        100
    FeedBin_k2_t1 OneTank_t1 1
    FeedBin_k2_t1 MaxRate_k1_t0 -120
    FeedBin_k2_t1 MinRate_k1_t0 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k2_t1 OBJ        14
    Flow_k2_t1 CDU_Demand_t1 1
    Flow_k2_t1 MaxRate_k1_t0 1
    Flow_k2_t1 MinRate_k1_t0 1
    Flow_k2_t1 TankInv_k2 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k2_t2 OBJ        100
    FeedBin_k2_t2 OneTank_t2 1
    FeedBin_k2_t2 MaxRate_k1_t1 -120
    FeedBin_k2_t2 MinRate_k1_t1 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k2_t2 OBJ        14
    Flow_k2_t2 CDU_Demand_t2 1
    Flow_k2_t2 MaxRate_k1_t1 1
    Flow_k2_t2 MinRate_k1_t1 1
    Flow_k2_t2 TankInv_k2 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k2_t3 OBJ        100
    FeedBin_k2_t3 OneTank_t3 1
    FeedBin_k2_t3 MaxRate_k1_t2 -120
    FeedBin_k2_t3 MinRate_k1_t2 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k2_t3 OBJ        14
    Flow_k2_t3 CDU_Demand_t3 1
    Flow_k2_t3 MaxRate_k1_t2 1
    Flow_k2_t3 MinRate_k1_t2 1
    Flow_k2_t3 TankInv_k2 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k2_t4 OBJ        100
    FeedBin_k2_t4 OneTank_t4 1
    FeedBin_k2_t4 MaxRate_k1_t3 -120
    FeedBin_k2_t4 MinRate_k1_t3 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k2_t4 OBJ        14
    Flow_k2_t4 CDU_Demand_t4 1
    Flow_k2_t4 MaxRate_k1_t3 1
    Flow_k2_t4 MinRate_k1_t3 1
    Flow_k2_t4 TankInv_k2 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k3_t1 OBJ        150
    FeedBin_k3_t1 OneTank_t1 1
    FeedBin_k3_t1 MaxRate_k2_t0 -120
    FeedBin_k3_t1 MinRate_k2_t0 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k3_t1 OBJ        18
    Flow_k3_t1 CDU_Demand_t1 1
    Flow_k3_t1 MaxRate_k2_t0 1
    Flow_k3_t1 MinRate_k2_t0 1
    Flow_k3_t1 TankInv_k3 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k3_t2 OBJ        150
    FeedBin_k3_t2 OneTank_t2 1
    FeedBin_k3_t2 MaxRate_k2_t1 -120
    FeedBin_k3_t2 MinRate_k2_t1 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k3_t2 OBJ        18
    Flow_k3_t2 CDU_Demand_t2 1
    Flow_k3_t2 MaxRate_k2_t1 1
    Flow_k3_t2 MinRate_k2_t1 1
    Flow_k3_t2 TankInv_k3 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k3_t3 OBJ        150
    FeedBin_k3_t3 OneTank_t3 1
    FeedBin_k3_t3 MaxRate_k2_t2 -120
    FeedBin_k3_t3 MinRate_k2_t2 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k3_t3 OBJ        18
    Flow_k3_t3 CDU_Demand_t3 1
    Flow_k3_t3 MaxRate_k2_t2 1
    Flow_k3_t3 MinRate_k2_t2 1
    Flow_k3_t3 TankInv_k3 1
    MARK0000  'MARKER'                 'INTORG'
    FeedBin_k3_t4 OBJ        150
    FeedBin_k3_t4 OneTank_t4 1
    FeedBin_k3_t4 MaxRate_k2_t3 -120
    FeedBin_k3_t4 MinRate_k2_t3 -80
    MARK0001  'MARKER'                 'INTEND'
    Flow_k3_t4 OBJ        18
    Flow_k3_t4 CDU_Demand_t4 1
    Flow_k3_t4 MaxRate_k2_t3 1
    Flow_k3_t4 MinRate_k2_t3 1
    Flow_k3_t4 TankInv_k3 1
RHS
    RHS1       OneTank_t1 1
    RHS1       CDU_Demand_t1 95
    RHS1       OneTank_t2 1
    RHS1       CDU_Demand_t2 95
    RHS1       OneTank_t3 1
    RHS1       CDU_Demand_t3 95
    RHS1       OneTank_t4 1
    RHS1       CDU_Demand_t4 95
    RHS1       TankInv_k1 200
    RHS1       TankInv_k2 180
    RHS1       TankInv_k3 150
RANGES
    RNG1       CDU_Demand_t1 25
    RNG1       CDU_Demand_t2 25
    RNG1       CDU_Demand_t3 25
    RNG1       CDU_Demand_t4 25
BOUNDS
 BV BND1       FeedBin_k1_t1
 UP BND1       Flow_k1_t1 120
 BV BND1       FeedBin_k1_t2
 UP BND1       Flow_k1_t2 120
 BV BND1       FeedBin_k1_t3
 UP BND1       Flow_k1_t3 120
 BV BND1       FeedBin_k1_t4
 UP BND1       Flow_k1_t4 120
 BV BND1       FeedBin_k2_t1
 UP BND1       Flow_k2_t1 120
 BV BND1       FeedBin_k2_t2
 UP BND1       Flow_k2_t2 120
 BV BND1       FeedBin_k2_t3
 UP BND1       Flow_k2_t3 120
 BV BND1       FeedBin_k2_t4
 UP BND1       Flow_k2_t4 120
 BV BND1       FeedBin_k3_t1
 UP BND1       Flow_k3_t1 120
 BV BND1       FeedBin_k3_t2
 UP BND1       Flow_k3_t2 120
 BV BND1       FeedBin_k3_t3
 UP BND1       Flow_k3_t3 120
 BV BND1       FeedBin_k3_t4
 UP BND1       Flow_k3_t4 120
ENDATA
