//
//  Pipeline.cpp
//
//  Created by Warren R. Carithers on 2016/10/19.
//  Based on a version created by Joe Geigel on 11/30/11.
//  Copyright 2016 Rochester Institute of Technology. All rights reserved.
//
//  Contributor:  Ruiyang Xu
//

#include <iostream>
#include "Pipeline.h"
#include <list>
using namespace std;
// Simple class that performs rasterization algorithms
//
// @param C The Canvas to use
///
Rasterizer::Rasterizer(Canvas &canvas ) : C(canvas)
{
}

///
// Get an edge table
//
// @param n - number of vertices
// @param x - array of x coordinates
// @param y - array of y coordinates
///
vector<list<Line> > Rasterizer::GetEdgeTable(int n, const int x[], const int y[])
{
    vector<list<Line> > edge_table;
    for (int i = 0; i <= Rasterizer::GetMax(n, y); i++) {
        list<Line> list;
        edge_table.push_back(list);
    }
    for (int i = 0; i < n; i++) {
            int j;
            if (i == n - 1) {
                j = 0;
            } else {
                j = i + 1;
            }
            if(y[j] != y[i]) {
                int ymin = min(y[i], y[j]);
                int ymax = max(y[i], y[j]);
                int min_index = i;
                if (y[i] > y[j]) {
                    min_index = j;
                }
                int max_index = j;
                if (y[i] > y[j]) {
                    max_index = i;
                }
                //double m_inverse = (y[j] == y[i] ? 0 : (double) (x[max_index] - x[min_index]) /
                                                       //(double) (y[max_index] - y[min_index]));
                double m_inverse = (double) (x[max_index] - x[min_index]) /(double) (y[max_index] - y[min_index]);
                Line line(x[min_index], ymax, m_inverse);
                edge_table[ymin].push_back(line);
            }
    }
    return edge_table;
}

///
// Get intersections from the sorted active list
//
// @param n - number of vertices
// @param y - array of y coordinates
///
vector<Range> get_ranges(list<Line>& active_list)
{
    vector<Range> ranges;
    if (!active_list.empty()) {
        list<Line>::iterator last_iter = active_list.end();
        advance(last_iter, -1);
        for (list<Line>::iterator iter = active_list.begin(); iter != last_iter; iter++) {
            list<Line>::iterator a_iterator = iter;
            advance(a_iterator, 1);
            Range r(iter->xmin, a_iterator->xmin);
            ranges.push_back(r);
        }
    }
    return ranges;
}


///
// Draw a filled polygon.
//
// Implementation should use the scan-line polygon fill algorithm
// discussed in class.
//
// The polygon has n distinct vertices.  The coordinates of the vertices
// making up the polygon are stored in the x and y arrays.  The ith
// vertex will have coordinate (x[i],y[i]).
//
// You are to add the implementation here using only calls to the
// setPixel() function.
//
// @param n - number of vertices
// @param x - array of x coordinates
// @param y - array of y coordinates
///
void Rasterizer::drawPolygon(int n, const int x[], const int y[] )
{
    list<Line> active_list;
    vector<list<Line > > edge_table = GetEdgeTable(n, x, y);

    //Drawing scanline by scanline
    for (int scan_line = GetMin(n, y); scan_line < edge_table.size(); scan_line++) {
        //get rid of unnecessary edge
        for (list<Line>::iterator iter = active_list.begin(); iter != active_list.end();) {
            if (iter->ymax <= scan_line) {
                iter = active_list.erase(iter);
            } else {
                iter++;
            }
        }
        active_list.splice(active_list.end(), edge_table[scan_line]);
        //sort active list on x intersection
        active_list.sort();
        // get x coordinate of ranges of intersections
        vector<Range> ranges = get_ranges(active_list);

        int xmin = GetMin(n, x);
        int xmax = GetMax(n, x);
        bool parity = true;
        for (unsigned int i = 0; i < ranges.size(); i++)
        {
            for (int j = ranges[i].min_range; j < ranges[i].max_range; j++)
            {
                if (parity) {
                    C.setPixel(j, scan_line);
                }
            }
            parity = !parity;
        }

        for (list<Line>::iterator iter = active_list.begin(); iter != active_list.end(); iter++)
        {
            iter->xmin += iter->inverse_m;
        }
    }
}

///
// Get min value from a list of values
//
// @param n - number of vertices
// @param y - array of y coordinates
///
int Rasterizer::GetMin(int n, const int y[])
{
    int min = y[0];
    for (int i = 0; i < n; i++) {
        if (y[i] < min) {
            min = y[i];
        }
    }
    return min;
}

// Get max value from a list of values
//
// @param n - number of vertices
// @param y - array of y coordinates
///
int Rasterizer::GetMax(int n, const int y[]) {
    int max = y[0];
    for (int i = 0; i < n; i++) {
        if (y[i] > max) {
            max = y[i];
        }
    }
    return max;
}

///
// Simple wrapper class for midterm assignment
//
// Only methods in the implementation file whose bodies
// contain the comment
//
//    // YOUR IMPLEMENTATION HERE
//
// are to be modified by students.
///


///
// Constructor
//
// @param w width of canvas
// @param h height of canvas
///
Pipeline::Pipeline( int w, int h ) : Canvas(w,h), mPolygons(), mClipper(), mViewPort(), mWorld(1.0f), mRasterizer(*this)
    // YOUR IMPLEMENTATION HERE if you need to add initializers
{
    // YOUR IMPLEMENTATION HERE if you need to modify the constructor
}

///
// addPoly - Add a polygon to the canvas.  This method does not draw
//           the polygon, but merely stores it for later drawing.
//           Drawing is initiated by the drawPoly() method.
//
//           Returns a unique integer id for the polygon.
//
// @param p - Array containing the vertices of the polygon to be added.
// @param n - Number of vertices in polygon
//
// @return a unique integer identifier for the polygon
///
int Pipeline::addPoly( const Vertex p[], int n )
{
    // YOUR IMPLEMENTATION HERE

    // REMEMBER TO RETURN A UNIQUE ID FOR THE POLYGON
	int id = mPolygons.size();
	vector<Vertex> poly;
	for (int i = 0; i < n; i++)
	{
		poly.push_back(p[i]);
	}
	mPolygons.insert(std::pair<int, vector<Vertex> >(id, poly));
    return id;
}

///
// drawPoly - Draw the polygon with the given id.  The polygon should
//            be drawn after applying the current transformation to
//            the vertices of the polygon.
//
// @param polyID - the ID of the polygon to be drawn.
///
void Pipeline::drawPoly( int polyID )
{
	vector<Vertex> poly = mPolygons[polyID];
	vector<Vertex> tPoly;
	for (unsigned int i = 0; i < poly.size(); i++)
	{
		Vertex vertex;
		glm::vec2 p = glm::vec2(mWorld * glm::vec3(poly[i].x, poly[i].y, 1.0f));
		vertex.x = p.x;
		vertex.y = p.y;
		tPoly.push_back(vertex);
	}
	// Clipping
	Vertex out[50];
	int n = mClipper.clipPolygon(tPoly.size(), tPoly.data(), out);


	//Viewport Transform
	float xd_min = mViewPort.xmin;
	float xd_max = mViewPort.xmin + mViewPort.width;
	float yd_min = mViewPort.ymin;
	float yd_max = mViewPort.ymin + mViewPort.height;

	glm::mat3 viewPortTransform(glm::vec3((xd_max - xd_min) / (mClipper.m_right - mClipper.m_left), 0.0f, 0.0f),
								glm::vec3(0.0f, (yd_max - yd_min) / (mClipper.m_top - mClipper.m_bottom), 0.0f),
								glm::vec3(
								(mClipper.m_right * xd_min - mClipper.m_left * xd_max) / (mClipper.m_right - mClipper.m_left),
								(mClipper.m_top * yd_min - mClipper.m_bottom * yd_max) / (mClipper.m_top - mClipper.m_bottom),
								1.0f)
								);

	int* x = new int[n];
	int* y = new int[n];
	for (unsigned int i = 0; i < n; i++)
	{
		glm::vec3 final_postition = viewPortTransform * glm::vec3(out[i].x, out[i].y, 1.0f);
		x[i] = final_postition.x;
		y[i] = final_postition.y;
	}

	// draw the final position
	mRasterizer.drawPolygon(n, x, y);

    //free memories
	delete[] x;
	delete[] y;
}

///
// clearTransform - Set the current transformation to the identity matrix.
///
void Pipeline::clearTransform( void )
{
	mWorld = glm::mat3(1.0f);
}

///
// translate - Add a translation to the current transformation by
//             premultiplying the appropriate translation matrix to
//             the current transformation matrix.
//
// @param x - Amount of translation in x.
// @param y - Amount of translation in y.
///
void Pipeline::translate( float tx, float ty )
{
	glm::mat3 trans(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(tx, ty, 1.0f));
	mWorld = trans * mWorld;
}

///
// rotate - Add a rotation to the current transformation by premultiplying
//          the appropriate rotation matrix to the current transformation
//          matrix.
//
// @param degrees - Amount of rotation in degrees.
///
void Pipeline::rotate( float degrees )
{
	glm::mat3 rotate(glm::vec3(
								glm::cos(glm::radians(degrees)), 
								glm::sin(glm::radians(degrees)),
								0.0f),
					 glm::vec3(
						 -1 * glm::sin(glm::radians(degrees)),
						 glm::cos(glm::radians(degrees)),
						 0.0f),
					glm::vec3(0.0f, 0.0f, 1.0f)
	);
	mWorld = rotate * mWorld;
}

///
// scale - Add a scale to the current transformation by premultiplying
//         the appropriate scaling matrix to the current transformation
//         matrix.
//
// @param x - Amount of scaling in x.
// @param y - Amount of scaling in y.
///
void Pipeline::scale( float sx, float sy )
{
    // YOUR IMPLEMENTATION HERE
	glm::mat3 scale(glm::vec3(sx, 0.0f,0.0f),
					glm::vec3(0.0f, sy, 0.0f),
					glm::vec3(0.0f, 0.0f, 1.0f)
	);
	mWorld = scale * mWorld;
}

///
// setClipWindow - Define the clip window.
//
// @param bottom - y coord of bottom edge of clip window (in world coords)
// @param top - y coord of top edge of clip window (in world coords)
// @param left - x coord of left edge of clip window (in world coords)
// @param right - x coord of right edge of clip window (in world coords)
///
void Pipeline::setClipWindow( float bottom, float top, float left, float right )
{
    // YOUR IMPLEMENTATION HERE
	mClipper.m_bottom = bottom;
	mClipper.m_top = top;
	mClipper.m_left = left;
	mClipper.m_right = right;
}

///
// setViewport - Define the viewport.
//
// @param xmin - x coord of lower left of view window (in screen coords)
// @param ymin - y coord of lower left of view window (in screen coords)
// @param width - width of view window (in world coords)
// @param height - width of view window (in world coords)
///
void Pipeline::setViewport( int x, int y, int width, int height )
{
    // YOUR IMPLEMENTATION HERE
	mViewPort.xmin = x;
	mViewPort.ymin = y;
	mViewPort.width = width - 1;
	mViewPort.height = height - 1;
}

///
// Clipper Constructor
///
Clipper::Clipper(float l, float r, float t, float b) 
	: m_left(l), m_right(r), m_top(t), m_bottom(b)
{
}

///
// Clipper Constructor
///
Clipper::Clipper()
	: m_left(0.0f), m_right(0.0f), m_top(0.0f), m_bottom(0.0f)
{
}

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
int Clipper::clipPolygon(int in, const Vertex inV[], Vertex outV[])
{
	// YOUR CODE GOES HERE
	Vertex out_left[50];
	Vertex out_right[50];
	Vertex out_top[50];
	int num_outputs = SHPC(in, inV, out_left, left, m_left);
	num_outputs = SHPC(num_outputs, out_left, out_right, right, m_right);
	num_outputs = SHPC(num_outputs, out_right, out_top, top, m_top);
	num_outputs = SHPC(num_outputs, out_top, outV, bottom, m_bottom);
	return num_outputs;  // remember to return the outgoing vertex count!
}

/*
 * Determine a vertex is either in or outside of the boundary
*/
bool Clipper::vertex_inside(Direction dir, float boundary, const Vertex & vertex)
{
	switch (dir)
	{
	case left:
		if (vertex.x < boundary)
		{
			return false;
		}
		else
		{
			return true;
		}
		break;
	case right:
		if (vertex.x > boundary)
		{
			return false;
		}
		else
		{
			return true;
		}
		break;
	case top:
		if (vertex.y > boundary)
		{
			return false;
		}
		else
		{
			return true;
		}
		break;
	case bottom:
		if (vertex.y < boundary)
		{
			return false;
		}
		else
		{
			return true;
		}
		break;
	}
}

/*
 * calculate the intersections of a boundary and a line segment given two vertices
 * @param dir    direction of a boundary edge
 * @param boundary   boundary
 * @param vertex1
 * @param vertex2
 * */
Vertex Clipper::intersection(Direction dir, float boundary, const Vertex & vertex1, const Vertex & vertex2)
{
	if (dir == left || dir == right)
	{
		float m = (vertex2.y - vertex1.y) / (vertex2.x - vertex1.x);
		float b = vertex1.y - m * vertex1.x;
		Vertex vertex;
		vertex.x = boundary;
		vertex.y = m * vertex.x + b;
		return vertex;
	}
	else
	{
		if (vertex2.x == vertex1.x)
		{
			Vertex vertex;
			vertex.y = boundary;
			vertex.x = vertex1.x;
			return vertex;
		}
		else
		{
			float m = (vertex2.y - vertex1.y) / (vertex2.x - vertex1.x);
			float b = vertex1.y - m * vertex1.x;
			Vertex vertex;
			vertex.y = boundary;
			vertex.x = (vertex.y - b) / m;
			return vertex;
		}
	}
}

/*
 * Sutherland Hodgman algorithm
 * @param number of input size
 * @param input vertices
 * @param direction of the boundary
 * @param boundary
 * */
int Clipper::SHPC(int n, const Vertex input[], Vertex output[], Direction dir, float boundary)
{
	int out_size = 0;
	const Vertex* end = &input[n - 1];
	for (int i = 0; i < n; i++)
	{
		const Vertex* start = &input[i];
		if (vertex_inside(dir, boundary, *start))
		{
			if (vertex_inside(dir, boundary, *end))
			{
				output[out_size++] = *start;
			}
			else
			{
				Vertex inter = intersection(dir, boundary, *start, *end);
				output[out_size++] = inter;
				output[out_size++] = *start;
			}
		}
		else
		{
			if (vertex_inside(dir, boundary, *end))
			{
				Vertex inter = intersection(dir, boundary, *start, *end);
				output[out_size++] = inter;
			}
		}
		end = start;
	}
	return out_size;
}
