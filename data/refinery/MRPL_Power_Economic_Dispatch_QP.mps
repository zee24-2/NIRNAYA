NAME          MRPL_Power_Economic_Dispatch_QP
OBJSENSE
  MIN
ROWS
 N  OBJ
 E  Power_Balance_180MW
COLUMNS
    Gen1_GT    OBJ        18
    Gen1_GT    Power_Balance_180MW 1
    Gen2_STG   OBJ        22
    Gen2_STG   Power_Balance_180MW 1
    Gen3_Aux   OBJ        26
    Gen3_Aux   Power_Balance_180MW 1
RHS
    RHS1       Power_Balance_180MW 180
RANGES
BOUNDS
 LO BND1       Gen1_GT    20
 UP BND1       Gen1_GT    90
 LO BND1       Gen2_STG   15
 UP BND1       Gen2_STG   80
 LO BND1       Gen3_Aux   10
 UP BND1       Gen3_Aux   70
QUADOBJ
    Gen1_GT    Gen1_GT    0.12
    Gen2_STG   Gen2_STG   0.16
    Gen3_Aux   Gen3_Aux   0.20000000000000001
ENDATA
