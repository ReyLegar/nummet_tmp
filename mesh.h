#pragma once

struct Point
{
    double x;
    double y;

    Point() : x(0.0), y(0.0) {}
    Point(double ax, double ay) : x(ax), y(ay) {}
    inline void operator *= (double q) { x *= q; y *= q; }
    inline void operator /= (double q) { x /= q; y /= q; }

    inline void operator  = (double q) { x = q;  y = q; }
    inline void operator  = (Point p) { x = p.x;  y = p.y; }

    inline void operator += (Point p) { x += p.x; y += p.y; }
    inline void operator -= (Point p) { x -= p.x; y -= p.y; }

    inline void operator += (double q) { x += q; y += q; }
    inline void operator -= (double q) { x -= q; y -= q; }

    friend inline double operator * (Point p1, Point p2) { return p1.x * p2.x + p1.y * p2.y; }
    friend inline Point operator + (Point p1, Point p2) { return Point(p1.x + p2.x, p1.y + p2.y); }
    friend inline Point operator - (Point p1, Point p2) { return Point(p1.x - p2.x, p1.y - p2.y); }
};

typedef Point Vector;


class Cell
{
public:
    Cell() : nodesInd(0), edgesInd(0), nCount(0), eCount(0) {};
    ~Cell();

    int   nCount;
    int   eCount;
    int* nodesInd;
    int* edgesInd;
    int	  neigh[3];
    int	  type;
    char  typeName[64];
    double  S;
    Point c; // центр ячейки
    double HX;
    double HY;
    unsigned int flag;

    friend class Grid;
};


class Edge
{
public:
    Edge() : c(0), cCount(0) {};
    ~Edge();

    int      n1;        // узел в начале
    int      n2;        // узел в конце
    int      c1;        // ячейка слева
    int      c2;        // ячейка справа, нормаль из с1 в с2
    Vector   n;         // нормаль к грани
    double     l;         // длина грани
    double     cnl1;       // сумма длин проекций векторов из центра грани до центра смежной ячейки
    double     cnl2;       // сумма длин проекций векторов из центра грани до центра смежной ячейки
    int      cCount;    // количество точек на грани
    Point* c;         // точки на грани
    int      type;      // тип грани (внутр., гранич.)
    friend class Grid;
public:
    static const int TYPE_INNER = 0;
    static const int TYPE_INLET = 1;
    static const int TYPE_OUTLET = 2;
    static const int TYPE_WALL = 3;

    static const int TYPE_NAMED = 0x100;
};



class Mesh
{
public:
    Mesh() : nodes(0), cells(0), edges(0),
        nCount(0), cCount(0), eCount(0) {};
    ~Mesh();

    void initFromFiles(char* fName);
    inline Point& getNode(int i) { return nodes[i]; };
    inline Cell& getCell(int i) { return cells[i]; };
    inline Edge& getEdge(int i) { return edges[i]; };
    int findEdge(int n1, int n2);

    Point* nodes;
    Cell* cells;
    Edge* edges;
    int    nCount;	//количество nodes.
    int    cCount;  //количество ячеек(cells).
    int    eCount;  //количество ребер(edges).
};