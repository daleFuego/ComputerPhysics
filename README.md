# ComputerPhysics

Zadanie 2. Pole elektryczne pomiędzy dwiema elektrodami

Oblicz pole elektryczne pomiędzy dwiema elektrodami – jedną kwadratową i jedną okrągłą. Wymiary elektrody (średnica a i długość boku b) oraz ich potencjały mogą być dobrane przez studenta. Elektroda okrągła znajduje się wewnątrz elektrody kwadratowej, centra obu nie muszą się pokrywać. Student powinien umieć zidentyfikować obszary o największych wartościach pola i przedstawić wykresy intensywności i potencjału. Dokładność rozwiązania musi być zweryfikowana. Rozwiązanie powinno zostać zaimplementowane dla minimum trzech rozmiarów siatki.

clear
xdel(winsid());

GRID_DENSITY = 8;
GRID_MIDDLE = (GRID_DENSITY / 2) + 1;

CIRCLE_ELECTRODE_RADIUS = 1;
CIRCLE_ELECTRODE_VALUE = 20;

CIRCLE_CENTRE_X = 0;
CIRCLE_CENTRE_Y = 0;

CONVERGENCE = 1e-6;
DIVERGENCE = 1e2;
MAX_CYCLES = 1e5;

GRID_AXIS_MIN = -CIRCLE_ELECTRODE_RADIUS;
GRID_AXIS_MAX = CIRCLE_ELECTRODE_RADIUS;

// 1. Tworzenie płaszczyzny z elektrodą okrągłą
moreIterations = %t;

x = linspace(GRID_AXIS_MIN, GRID_AXIS_MAX, GRID_DENSITY);
y = x;
gridStep = x(2)-x(1);

[X,Y] = ndgrid(x,x);
circleMatrix = (X - CIRCLE_CENTRE_X).^2 + (Y - CIRCLE_CENTRE_Y).^2 < CIRCLE_ELECTRODE_RADIUS^2;

// Szukamy pól w jakich możemy liczyć pole ?
LN = [circleMatrix(:,$),circleMatrix(:,1:$-1)];
RN = [circleMatrix(:,2:$),circleMatrix(:,1)];
TN = [circleMatrix($,:);circleMatrix(1:$-1,:)];
BN = [circleMatrix(2:$,:);circleMatrix(1,:)];

// Dostajemy mapę zasięgu elektrody ?
numberOfNeighbours = LN+RN+TN+BN;

// Pomocnicza macierz zawierająca same false'y, bo wszędzie w X mamy wartości liczbowe, a nie Not A Number
nAnDetection = isnan(X);

// Szukamy krawędzi prostych
// single top nneighbour
top = nAnDetection;
top(:, GRID_MIDDLE:$-1) = (numberOfNeighbours(:,GRID_MIDDLE:$-1)==3) & ~circleMatrix(:,GRID_MIDDLE+1:$);
// single bottom neighbour
bottom = nAnDetection;
bottom(:,2:GRID_MIDDLE+1) = (numberOfNeighbours(:,2:GRID_MIDDLE+1)==3) & ~circleMatrix(:, 1:GRID_MIDDLE);
// single left neighbour
left = nAnDetection;
left(2:GRID_MIDDLE+1, :) = (numberOfNeighbours(2:GRID_MIDDLE+1,:)==3) & ~circleMatrix(1:GRID_MIDDLE,:);
// single right neighbour
right = nAnDetection;
right(GRID_MIDDLE:$-1,:) = (numberOfNeighbours(GRID_MIDDLE:$-1,:)==3) & ~circleMatrix(GRID_MIDDLE+1:$,:);

// Exclude TOP RADIUS point (middle on)!

// Szukamy krawędzi zagiętych
bottomLeft = nAnDetection;
bottomLeft(:,2:GRID_MIDDLE-1) = (numberOfNeighbours(:,2:GRID_MIDDLE-1)==2) & circleMatrix(:,2:GRID_MIDDLE-1);
bottomLeft(:,2:GRID_MIDDLE-1) = bottomLeft(:,2:GRID_MIDDLE-1) & bottomLeft(2:GRID_MIDDLE-1,:)';

bottomRight = nAnDetection;
bottomRight(:,2:GRID_MIDDLE-1) = (numberOfNeighbours(:,2:GRID_MIDDLE-1)==2) & circleMatrix(:,2:GRID_MIDDLE-1);
// Pertranspose ? :D ?
temp = pertrans(bottomRight);
bottomRight = bottomRight & temp;

topLeft = nAnDetection;
topLeft(:, GRID_MIDDLE+1:$-1) = (numberOfNeighbours(:,GRID_MIDDLE+1:$-1)==2) & circleMatrix(:,GRID_MIDDLE+1:$-1);
temp = pertrans(topLeft);
topLeft(:, GRID_MIDDLE+1:$-1) = topLeft(:, GRID_MIDDLE+1:$-1) & temp(:, GRID_MIDDLE+1:$-1);

topRight = nAnDetection;
topRight(:, GRID_MIDDLE+1:$-1) = (numberOfNeighbours(:,GRID_MIDDLE+1:$-1)==2) & circleMatrix(:,GRID_MIDDLE+1:$-1);
topRight = topRight & topRight';

// Finding the special points:
// boundary on the right
[xR,yR] = find(right);
nR = length(xR);
// boundary on the left
[xL,yL] = find(left);
nL = length(xL);
// boundary on the top
[xT,yT] = find(top);
nT = length(xT);
// boundary on the bottom
[xB,yB] = find(bottom);
nB = length(xB);
// boundary on the top-right
[xTR,yTR] = find(topRight);
nTR = length(xTR);
// boundary on the top-left
[xTL,yTL] = find(topLeft);
nTL = length(xTL);
// boundary on the bottom-right
[xBR,yBR] = find(bottomRight);
nBR = length(xBR);
// boundary on the bottom-left
[xBL,yBL] = find(bottomLeft);
nBL = length(xBL);

// Vectors with boundary coordinates:
xb = sqrt(CIRCLE_ELECTRODE_RADIUS^2 - (y).^2);
yb = sqrt(CIRCLE_ELECTRODE_RADIUS^2 - (x).^2);

// Determination of distances of special points from the boundary:
nr_xb = sum(X<ones(x')*xb,"c")';
delta_xb = [xb(1:$-1) - x(nr_xb(1:$-1)) , 0];
nr_yb = sum(Y<ones(y)*yb','c')';
delta_yb = [yb(1:$-1) - y(nr_yb(1:$-1)) , 0];

// coefficients "a" and "b":
a = delta_xb / gridStep;
temp1 = pertrans(a)';
a(1:$/2)= temp1(1:$/2);

b = delta_yb / gridStep;
b =-b;
temp2 = pertrans(b)';
b(1:$/2)= temp2(1:$/2);

function [initialField]=createInitialField(initialField)
    circumference = circleMatrix;
    initialField (circumference == %T) = CIRCLE_ELECTRODE_VALUE;
endfunction

ininitialField = zeros(X)
ininitialField = createInitialField(ininitialField);
newField = ininitialField;

while moreIterations do
    newField(2:$-1,2:$-1) = [ininitialField(1:$-2,2:$-1) + ininitialField(3:$,2:$-1) + ininitialField(2:$-1,1:$-2) + ininitialField(2:$-1,3:$)]/4;

    for i = 1:nR do
        newField(xR(i),yR(i)) = (CIRCLE_ELECTRODE_VALUE+a(xR(i))*ininitialField(xR(i)-1,yR(i))) /..
        (a(xR(i))+1);
    end
    
    // left
    for i = 1:nL do
        newField(xL(i),yL(i)) = (CIRCLE_ELECTRODE_VALUE+a(xL(i))*ininitialField(xL(i)+1,yL(i))) /..
        (a(xL(i))+1);
    end
    
    // top
    for i = 1:nT do
        newField(xT(i),yT(i)) = (CIRCLE_ELECTRODE_VALUE+b(yT(i))*ininitialField(xT(i),yT(i)-1)) /..
        (b(yT(i))+1);
    end
    
    // bottom
    for i = 1:nB do
        newField(xB(i),yB(i)) = (CIRCLE_ELECTRODE_VALUE+b(yB(i))*ininitialField(xB(i),yB(i)+1)) /..
        (b(yB(i))+1);
    end
    
    //top-right
    for i = 1:nTR do
        numerator = (ininitialField(xTR(i)-1,yTR(i))*a(xTR(i))*b(yTR(i))*(b(yTR(i))+1))+..
        (CIRCLE_ELECTRODE_VALUE*b(yTR(i))*(b(yTR(i))+1))+(CIRCLE_ELECTRODE_VALUE*a(xTR(i))*(a(xTR(i))+1))+..
        (ininitialField(xTR(i),yTR(i)-1)*a(xTR(i))*b(yTR(i))*(a(xTR(i))+1));
        denominator = (a(xTR(i))+1) * (b(yTR(i))+1) * (a(xTR(i)) + b(yTR(i)));
        newField(xTR(i),yTR(i)) = numerator / denominator;
    end
    
    //top-left
    for i = 1:nTL do
        numerator = (ininitialField(xTL(i)+1,yTL(i))*a(xTL(i))*b(yTL(i))*(b(yTL(i))+1))+..
        (CIRCLE_ELECTRODE_VALUE*b(yTL(i))*(b(yTL(i))+1))+(CIRCLE_ELECTRODE_VALUE*a(xTL(i))*(a(xTL(i))+1))+..
        (ininitialField(xTL(i),yTL(i)-1)*a(xTL(i))*b(yTL(i))*(a(xTL(i))+1));
        denominator = (a(xTL(i))+1) * (b(yTL(i))+1) * (a(xTL(i)) + b(yTL(i)));
        newField(xTL(i),yTL(i)) = numerator / denominator;
    end
    
    //bottom-right
    for i = 1:nBR do
        numerator = (ininitialField(xBR(i)-1,yBR(i))*a(xBR(i))*b(yBR(i))*(b(yBR(i))+1))+..
        (CIRCLE_ELECTRODE_VALUE*b(yBR(i))*(b(yBR(i))+1))+(CIRCLE_ELECTRODE_VALUE*a(xBR(i))*(a(xBR(i))+1))+..
        (ininitialField(xBR(i),yBR(i)+1)*a(xBR(i))*b(yBR(i))*(a(xBR(i))+1));
        denominator = (a(xBR(i))+1) * (b(yBR(i))+1) * (a(xBR(i)) + b(yBR(i)));
        newField(xBR(i),yBR(i)) = numerator / denominator;
    end
    
    //bottom-left
    for i = 1:nBL do
        numerator = (ininitialField(xBL(i)+1,yBL(i))*a(xBL(i))*b(yBL(i))*(b(yBL(i))+1))+..
        (CIRCLE_ELECTRODE_VALUE*b(yBL(i))*(b(yBL(i))+1))+(CIRCLE_ELECTRODE_VALUE*a(xBL(i))*(a(xBL(i))+1))+..
        (ininitialField(xBL(i),yBL(i)+1)*a(xBL(i))*b(yBL(i))*(a(xBL(i))+1));
        denominator = (a(xBL(i))+1) * (b(yBL(i))+1) * (a(xBL(i)) + b(yBL(i)));
        newField(xBL(i),yBL(i)) = numerator / denominator;
    end
    
    finish = max(abs(ininitialField-newField));
    if (finish < CONVERGENCE | finish > DIVERGENCE) then
         moreIterations = %f;
    end
    ininitialField = newField;
end

// 2. Tworzenie płaszczyzny z elektrodą kradratową
SQ_GRID_DENSITY = 4 * GRID_DENSITY;
SQ_GRID_MIDDLE = (SQ_GRID_DENSITY / 2) + 1;

SQUARE_ELECTRODE_SIZE = 16;
SQUARE_ELECTRODE_VALUE = -30;
SQUARE_START_X = SQ_GRID_MIDDLE - 5;
SQUARE_START_Y = SQ_GRID_MIDDLE - 5;

SQ_GRID_AXIS_MIN = -SQUARE_ELECTRODE_SIZE;
SQ_GRID_AXIS_MAX = SQUARE_ELECTRODE_SIZE;

CIRCLE_START_X = SQUARE_START_X + 7;
CIRCLE_START_Y = SQUARE_START_Y + 7;

x = linspace(SQ_GRID_AXIS_MIN, SQ_GRID_AXIS_MAX, SQ_GRID_DENSITY);
y = x;
z = zeros(SQ_GRID_DENSITY,SQ_GRID_DENSITY);

for idx = 1:SQUARE_ELECTRODE_SIZE
    for idy = 1:SQUARE_ELECTRODE_SIZE
        z(SQUARE_START_X + idx, SQUARE_START_Y + idy) = SQUARE_ELECTRODE_VALUE;
    end
end

for idx = 1:GRID_DENSITY
    for idy = 1:GRID_DENSITY
        z(CIRCLE_START_X + idx, CIRCLE_START_Y + idy) = z(CIRCLE_START_X + idx, CIRCLE_START_Y + idy) + newField(idx,idy);
    end
end

xset('colormap',jetcolormap(64));
plot3d1(x,y,z,alpha=5.5,theta=22.5);
colorbar(min(SQUARE_ELECTRODE_VALUE),max(CIRCLE_ELECTRODE_VALUE));
tyt = msprintf("Square electrode with circle electrode inside");
xtitle(tyt);
