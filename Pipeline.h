//
//  Pipeline.h
//
//  Created by Warren R. Carithers on 2016/10/19.
//  Based on a version created by Joe Geigel on 11/30/11.
//  Copyright 2016 Rochester Institute of Technology. All rights reserved.
//
//  Contributor:  Ruiyang Xu
//

#ifndef _PIPELINE_H_
#define _PIPELINE_H_

#include "Canvas.h"
#include "Vertex.h"
#include <map>
#include <list>
#include <glm/glm.hpp>
using namespace std;

class Canvas;

/// line bucket
struct Line
{
    double xmin;
    int ymax;
    double inverse_m;
    Line(double x, int y, double inverse):
            xmin(x), ymax(y), inverse_m(inverse)
    {
    }
    bool operator< (const Line& other)
    {
        return xmin < other.xmin;
    }
};

/// Intersection Range
struct Range
{
    double min_range;
    double max_range;
    Range(double m, double n) :min_range(m), max_range(n)
    {
    }

    bool is_in_range(int x) {
        return x > min_range && x < max_range;
    }
};


///
// Simple class that performs rasterization algorithms
///

class Rasterizer {
public:

    ///
    // Drawing canvas
    ///

    Canvas &C;

    ///
    // Constructor
    //
    // @param n number of scanlines
    // @param C The Canvas to use
    ///
    Rasterizer(Canvas &canvas );

    ///
    // Draw a filled polygon
    //
    // Implementation should use the scan-line polygon fill algorithm
    // discussed in class.
    //
    // The polygon has n distinct vertices.  The coordinates of the vertices
    // making up the polygon are stored in the x and y arrays.  The ith
    // vertex will have coordinate (x[i],y[i]).
    //
    // You are to add the implementation here using only calls to the
    // setPixel() method of the canvas.
    //
    // @param n - number of vertices
    // @param x - array of x coordinates
    // @param y - array of y coordinates
    ///
    void drawPolygon( int n, const int x[], const int y[] );
    // Get minimum value from a list of values
    //
    // @param n - number of vertices
    // @param y - array of y coordinates
    /// \param n
    /// \param y
    /// \return

private:
	static int GetMin(int n, const int y[]);

    // Get max value from a list of values
    //
    // @param n - number of vertices
    // @param y - array of y coordinates
    ///
    static int GetMax(int n, const int y[]);
    ///
    // Get an edge table
    //
    // @param n - number of vertices
    // @param x - array of x coordinates
    // @param y - array of y coordinates
    ///
    static vector<std::list<Line> > GetEdgeTable(int n, const int x[], const int y[]);
};

///
// Simple module that performs clipping
///
class Clipper
{
public:
    /// Direction
	enum Direction
	{
		left,
		right,
		top,
		bottom
	};
    ///
    // Constructor
    ///
	Clipper(float l, float r, float t, float b);
    ///
    // Constructor
    ///
	Clipper();
    /// Pipeline is the boss, it can access everything in Clipper
	friend class Pipeline;
    ///
    // clipPolygon
    //
    // Clip the polygon with vertex count in and vertices inV against the
    // rectangular clipping region specified by lower-left corner ll and
    // upper-right corner ur. The resulting vertices are placed in outV.
    //
    // The routine should return the with the vertex count of polygon
    // resulting from the clipping.
    //
    // @param in    the number of vertices in the polygon to be clipped
    // @param inV   the incoming vertex list
    // @param outV  the outgoing vertex list
    // @param ll    the lower-left corner of the clipping rectangle
    // @param ur    the upper-right corner of the clipping rectangle
    //
    // @return number of vertices in the polygon resulting after clipping
    //
    ///
	int clipPolygon(int in, const Vertex inV[], Vertex outV[]);
private:
	float m_left;
	float m_right;
	float m_top;
	float m_bottom;
	bool vertex_inside(Direction dir, float boundary, const Vertex& vertex);
	Vertex intersection(Direction dir, float boundary, const Vertex& vertex1, const Vertex& vertex2);
	int SHPC(int n, const Vertex input[], Vertex output[], Direction dir, float boundary);
};

///View port Description
struct ViewPort
{
	int xmin;
	int ymin;
	int width;
	int height;
};

class Pipeline : public Canvas {

public:

    ///
    // Constructor
    //
    // @param w width of canvas
    // @param h height of canvas
    ///
    Pipeline(int w, int h);

    ///
    // addPoly - Add a polygon to the canvas.  This method does not draw
    //           the polygon, but merely stores it for later drawing.
    //           Drawing is initiated by the drawPoly() method.
    //
    // @param p - Array of Vertex entries defining the polygon to be added
    // @param n - Number of vertices in polygon
    //
    // @return a unique integer identifier for the polygon
    ///
    int addPoly( const Vertex p[], int n );

    ///
    // drawPoly - Draw the polygon with the given id.  The polygon should
    //            be drawn after applying the current transformation to
    //            the vertices of the polygon.
    //
    // @param polyID - the ID of the polygon to be drawn.
    ///
    void drawPoly( int polyID );

    ///
    //
    // clearTransform - Set the current transformation to the identity matrix.
    //
    ///
    void clearTransform( void );

    ///
    // translate - Add a translation to the current transformation by
    //             premultiplying the appropriate translation matrix to
    //             the current transformation matrix.
    //
    // @param tx - Amount of translation in x.
    // @param ty - Amount of translation in y.
    //
    ///
    void translate( float tx, float ty );

    ///
    // rotate - Add a rotation to the current transformation by premultiplying
    //          the appropriate rotation matrix to the current transformation
    //          matrix.
    //
    // @param degrees - Amount of rotation in degrees.
    ///
    void rotate( float degrees );

    ///
    // scale - Add a scale to the current transformation by premultiplying
    //         the appropriate scaling matrix to the current transformation
    //         matrix.
    //
    // @param sx - Amount of scaling in x.
    // @param sy - Amount of scaling in y.
    ///
    void scale( float sx, float sy );

    ///
    // setClipWindow - Define the clip window.
    //
    // @param bottom - y coord of bottom edge of clip window (in world coords)
    // @param top - y coord of top edge of clip window (in world coords)
    // @param left - x coord of left edge of clip window (in world coords)
    // @param right - x coord of right edge of clip window (in world coords)
    ///
    void setClipWindow( float bottom, float top, float left, float right );

    ///
    // setViewport - Define the viewport.
    //
    // @param xmin - x coord of lower left of view window (in screen coords)
    // @param ymin - y coord of lower left of view window (in screen coords)
    // @param width - width of view window (in world coords)
    // @param height - width of view window (in world coords)
    ///
    void setViewport( int x, int y, int width, int height );
private:
    /// a clipper
	Clipper mClipper;
    /// polygons
	map<int, vector<Vertex> > mPolygons;
    /// view port description
	ViewPort mViewPort;
    /// model->world matrix
    glm::mat3 mWorld;
    /// a Rasterizer
    Rasterizer mRasterizer;
};

#endif
