

ans =

    0.0559

prop

prop =

    0.0558    0.0653    0.0599    0.0531    0.0563    0.0552    0.0535    0.0486

hist(max(CellR))
max(CellR)

ans =

    0.9583    0.8293    0.8438    0.7843    0.8333    0.7241    0.8636    0.7368

cell = find(CellR==0.9583)

cell =

  0×1 empty double column vector

CellR2= maxCellR([],2)
Undefined function or variable 'maxCellR'.
 
CellR2 = max(transpose(CellR))

CellR2 =

  Columns 1 through 11

    0.3158    0.3182    0.6897    0.3448    0.3902    0.3333    0.3333    0.5833    0.7059    0.3902    0.7917

  Columns 12 through 22

    0.2083    0.1667    0.0588    0.8293    0.2683    0.6829    0.0833    0.8049    0.1951    0.7843    0.5000

  Columns 23 through 33

    0.6078    0.1579    0.1379    0.4483    0.4167    0.2759    0.8125    0.4118    0.3448    0.7241    0.2500

  Columns 34 through 44

    0.2353    0.7188    0.5000    0.7805    0.6341    0.2745    0.4510    0.3137    0.0909    0.3171    0.2500

  Columns 45 through 55

    0.0980    0.0980    0.0784    0.4483    0.1951    0.2727    0.6562    0.1765    0.4375    0.6275    0.6562

  Columns 56 through 66

    0.2812    0.7500    0.5000    0.4211    0.7500    0.8333    0.6552    0.2439    0.5366    0.1569    0.0690

  Columns 67 through 77

    0.2745    0.3333    0.5490    0.7255    0.4091    0.2188    0.6207         0    0.3684    0.1667    0.4737

  Columns 78 through 88

    0.4211    0.2632    0.3158    0.1579    0.2105    0.0196    0.8182    0.8636    0.5455    0.4545    0.2500

  Columns 89 through 99

    0.7368    0.8333         0    0.6250    0.7500    0.7083    0.4167    0.7500    0.2105    0.0909    0.1250

  Columns 100 through 110

    0.3750    0.5909    0.4167    0.6250    0.6667    0.9583    0.7500    0.6667    0.4545    0.5000    0.5455

  Columns 111 through 121

    0.2500    0.0345    0.0833    0.0244         0    0.1373         0    0.8438    0.3684    0.2727    0.5789

  Columns 122 through 132

    0.1818    0.7368    0.5789    0.6863    0.6078    0.5833    0.2745    0.1562    0.3438    0.6316    0.2105

  Columns 133 through 143

    0.3158    0.3684    0.7500    0.5000    0.8333    0.2759    0.2439         0    0.0909    0.1707    0.1053

  Columns 144 through 149

    0.1579    0.0526    0.0417    0.0244    0.0833    0.0833

max (CellR2)

ans =

    0.9583

max(CellR2, index)
Undefined function or variable 'index'.
 
max = ([i],CellR2)
 max = ([i],CellR2)
           ?
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
max[m,i]=max(CellR2)
 max[m,i]=max(CellR2)
    ?
Error: Unbalanced or unexpected parenthesis or bracket.
 
[m,i]=max(CellR2)

m =

    0.9583


i =

   105

load('Cells.mat')
imagesc(XCellsCl)
Error using image
Image CData must be an m-by-n-by-3 or m-by-n matrix.

Error in imagesc (line 18)
    hh = image(varargin{1},'CDataMapping','scaled');
 
imagesc(XCellsCl,:,:,104)
Undefined function or variable 'imagesc'.
 
imagesc(XCellsCl(:,:,104))
imagesc(XCellsCl(:,:,105))
