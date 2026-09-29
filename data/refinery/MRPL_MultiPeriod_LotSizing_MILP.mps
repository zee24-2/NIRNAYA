NAME          MRPL_MultiPeriod_LotSizing_MILP
OBJSENSE
  MIN
ROWS
 N  OBJ
 E  Bal_t1
 L  CapLink_t1
 E  Bal_t2
 L  CapLink_t2
 E  Bal_t3
 L  CapLink_t3
 E  Bal_t4
 L  CapLink_t4
COLUMNS
    Prod_t1    OBJ        12
    Prod_t1    Bal_t1     1
    Prod_t1    CapLink_t1 1
    Inv_t1     OBJ        2.5
    Inv_t1     Bal_t1     -1
    Inv_t1     Bal_t2     1
    MARK0000  'MARKER'                 'INTORG'
    Setup_t1   OBJ        450
    Setup_t1   CapLink_t1 -250
    MARK0001  'MARKER'                 'INTEND'
    Prod_t2    OBJ        12
    Prod_t2    Bal_t2     1
    Prod_t2    CapLink_t2 1
    Inv_t2     OBJ        2.5
    Inv_t2     Bal_t2     -1
    Inv_t2     Bal_t3     1
    MARK0000  'MARKER'                 'INTORG'
    Setup_t2   OBJ        450
    Setup_t2   CapLink_t2 -250
    MARK0001  'MARKER'                 'INTEND'
    Prod_t3    OBJ        12
    Prod_t3    Bal_t3     1
    Prod_t3    CapLink_t3 1
    Inv_t3     OBJ        2.5
    Inv_t3     Bal_t3     -1
    Inv_t3     Bal_t4     1
    MARK0000  'MARKER'                 'INTORG'
    Setup_t3   OBJ        450
    Setup_t3   CapLink_t3 -250
    MARK0001  'MARKER'                 'INTEND'
    Prod_t4    OBJ        12
    Prod_t4    Bal_t4     1
    Prod_t4    CapLink_t4 1
    Inv_t4     OBJ        2.5
    Inv_t4     Bal_t4     -1
    MARK0000  'MARKER'                 'INTORG'
    Setup_t4   OBJ        450
    Setup_t4   CapLink_t4 -250
    MARK0001  'MARKER'                 'INTEND'
RHS
    RHS1       Bal_t1     120
    RHS1       Bal_t2     180
    RHS1       Bal_t3     150
    RHS1       Bal_t4     210
RANGES
BOUNDS
 UP BND1       Prod_t1    250
 UP BND1       Inv_t1     300
 BV BND1       Setup_t1
 UP BND1       Prod_t2    250
 UP BND1       Inv_t2     300
 BV BND1       Setup_t2
 UP BND1       Prod_t3    250
 UP BND1       Inv_t3     300
 BV BND1       Setup_t3
 UP BND1       Prod_t4    250
 UP BND1       Inv_t4     300
 BV BND1       Setup_t4
ENDATA
