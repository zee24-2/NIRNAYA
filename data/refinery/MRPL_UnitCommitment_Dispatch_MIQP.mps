NAME          MRPL_UnitCommitment_Dispatch_MIQP
OBJSENSE
  MIN
ROWS
 N  OBJ
 E  Load_100MW
 L  Cap_GT1
 L  Cap_GT2
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    Commit_GT1 OBJ        200
    Commit_GT1 Cap_GT1    -80
    Commit_GT2 OBJ        150
    Commit_GT2 Cap_GT2    -80
    MARK0001  'MARKER'                 'INTEND'
    Gen_GT1    OBJ        15
    Gen_GT1    Load_100MW 1
    Gen_GT1    Cap_GT1    1
    Gen_GT2    OBJ        18
    Gen_GT2    Load_100MW 1
    Gen_GT2    Cap_GT2    1
RHS
    RHS1       Load_100MW 100
RANGES
BOUNDS
 BV BND1       Commit_GT1
 BV BND1       Commit_GT2
 UP BND1       Gen_GT1    80
 UP BND1       Gen_GT2    80
QUADOBJ
    Gen_GT1    Gen_GT1    0.10000000000000001
    Gen_GT2    Gen_GT2    0.12
ENDATA
