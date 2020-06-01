#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    // Base Case - We are at the evaluated Bezier Curve point
      int numpts = points.size();
      std::vector<Vector2D> LERPpoints;
      if (numpts == 1) {
          return points;
      }
      // LERP-ing the list of points
      for (int i = 0; i < numpts-1; i++) {
          LERPpoints.push_back((1.0 - t) * points[i] + t * points[i + 1]);
          }
      return LERPpoints;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      // Base Case - We are at the evaluated Bezier Curve point
      int numpts = points.size();
      std::vector<Vector3D> LERPpoints;
      if (numpts == 1) {
          return points;
      }
      // LERP-ing the list of points
      for (int i = 0; i < numpts - 1; i++) {
          LERPpoints.push_back((1.0 - t) * points[i] + t * points[i + 1]);
      }
      return LERPpoints;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      // Base Case
      if (points.size() == 1) {
          return points[0];
      }
      return evaluate1D(evaluateStep(points, t), t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
      std::vector<Vector3D> movingBezierControlPoints;
      int n = controlPoints.size(); // outer dimension
      for (int i = 0; i < n ; i++) {
      movingBezierControlPoints.push_back(evaluate1D(controlPoints[i], u));      
      }

      return evaluate1D(movingBezierControlPoints, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.


      Vector3D weighted_normal(0,0,0);
      HalfedgeCIter h = halfedge();

      do {
              // don't add the area of boundary faces!
              if (!h->face()->isBoundary())
               {
                  // Get the area of the face
                  Vector3D v0 = position;
                  Vector3D v1 = h->next()->vertex()->position;
                  Vector3D v2 = h->next()->next()->vertex()->position;

                  double face_area = (cross(v1-v0,v2-v0).norm())/2;
                  weighted_normal += face_area * h->face()->normal();
          }
          h = h->twin()->next();
          
      }
      while (h != halfedge());

    return weighted_normal.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
      
      if (!e0->isBoundary()){
      
      // ---------------- DEFINITIONS ---------------------

      // Inner Halfedges
      HalfedgeIter h0 = e0->halfedge(); // halfedge 1
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();

      // Outer Halfedges
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      // Vertices
      VertexIter v0 = h0->vertex();  // first vertex at endpoint of previous edge
      VertexIter v1 = h3->vertex();  // second vertex at opposite endpoint of previous edge
      VertexIter v2 = h2->vertex();  // first vertex at endpoint of flipped edge
      VertexIter v3 = h5->vertex();  // second vertex at opposite endpoint of flipped edge

      // Edges
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();

      // Faces
      FaceIter f1 = h0->face();       // face 1
      FaceIter f2 = h3->face();       // face 2


      // --------------- UPDATES ------------------
      
      //Copy of set neighbors for argument list and order
      /*  void setNeighbors( 
          HalfedgeIter next,
          HalfedgeIter twin,
          VertexIter vertex,
          EdgeIter edge,
          FaceIter face ) 
      */

      // Inner Halfedges
      h0->setNeighbors(h1,h3,v2,e0,f1);
      h1->setNeighbors(h2,h9,v3,e4,f1);
      h2->setNeighbors(h0,h6,v1,e1,f1);
      h3->setNeighbors(h4,h0,v3,e0,f2);
      h4->setNeighbors(h5,h7,v2,e2,f2);
      h5->setNeighbors(h3,h8,v0,e3,f2);

      // Outer Halfedges
      h6->setNeighbors(h6->next(),h2,v2,e1,h6->face());
      h7->setNeighbors(h7->next(),h4,v0,e2,h7->face());
      h8->setNeighbors(h8->next(),h5,v3,e3,h8->face());
      h9->setNeighbors(h9->next(),h1,v1,e4,h9->face());

      // Vertices
      v0->halfedge() = h5;
      v1->halfedge() = h2;
      v2->halfedge() = h4;
      v3->halfedge() = h1;
      
      // Edges
      e0->halfedge() = h0;
      e1->halfedge() = h2;
      e2->halfedge() = h4;
      e3->halfedge() = h5;
      e4->halfedge() = h1;

      // Faces
      f1->halfedge() = h0;
      f2->halfedge() = h3;

      }
 

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
      if (!e0->isBoundary()) {
      // ---------------- DEFINITIONS ---------------------

      // Inner Halfedges
          HalfedgeIter h0 = e0->halfedge();
          HalfedgeIter h1 = h0->next();
          HalfedgeIter h2 = h1->next();
          HalfedgeIter h3 = h0->twin();
          HalfedgeIter h4 = h3->next();
          HalfedgeIter h5 = h4->next();

          // Outer Halfedges
          HalfedgeIter h6 = h1->twin();
          HalfedgeIter h7 = h2->twin();
          HalfedgeIter h8 = h4->twin();
          HalfedgeIter h9 = h5->twin();

          // Vertices
          VertexIter v0 = h0->vertex();  // first vertex at endpoint of previous edge
          VertexIter v1 = h3->vertex();  // second vertex at opposite endpoint of previous edge
          VertexIter v2 = h2->vertex();  // first vertex at endpoint of flipped edge
          VertexIter v3 = h5->vertex();  // second vertex at opposite endpoint of flipped edge

          // Edges
          EdgeIter e1 = h1->edge();
          EdgeIter e2 = h2->edge();
          EdgeIter e3 = h4->edge();
          EdgeIter e4 = h5->edge();

          // Faces
          FaceIter f1 = h0->face();       // face 1
          FaceIter f2 = h3->face();       // face 2

          // --------------- NEW ELEMENTS ------------
          // New inner halfedges
          HalfedgeIter h10 = newHalfedge();
          HalfedgeIter h11 = newHalfedge();
          HalfedgeIter h12 = newHalfedge();
          HalfedgeIter h13 = newHalfedge();
          HalfedgeIter h14 = newHalfedge();
          HalfedgeIter h15 = newHalfedge();

          // New vertex
          VertexIter v = newVertex();
          
          // New edges
          EdgeIter e5 = newEdge();
          EdgeIter e6 = newEdge();
          EdgeIter e7 = newEdge();

          // New faces
          FaceIter f3 = newFace();
          FaceIter f4 = newFace();


          // --------------- UPDATES ------------------

          // Halfedges
          h0->setNeighbors(h1,h3,v,e0,f1);
          h1->setNeighbors(h2,h6,v1,e1,f1);
          h2->setNeighbors(h0,h11,v2,e5,f1);
          h3->setNeighbors(h4,h0,v1,e0,f2);
          h4->setNeighbors(h5,h15,v,e7,f2);
          h5->setNeighbors(h3,h9,v3,e4,f2);
          h6->setNeighbors(h6->next(),h1,v2,e1,h6->face());
          h7->setNeighbors(h7->next(),h12,v0,e2,h7->face());
          h8->setNeighbors(h8->next(),h14,v3,e3,h8->face());
          h9->setNeighbors(h9->next(),h5,v1,e4,h9->face());
          h10->setNeighbors(h11,h13,v0,e6,f3);
          h11->setNeighbors(h12,h2,v,e5,f3);
          h12->setNeighbors(h10,h7,v2,e2,f3);
          h13->setNeighbors(h14,h10,v,e6,f4);
          h14->setNeighbors(h15,h8,v0,e3,f4);
          h15->setNeighbors(h13,h4,v3,e7,f4);

          // Vertices
          // Update vertex position to edge midpoint
          v->position = 0.5 * (v0->position + v1->position);
          // Update vertex status to isNew
          v->isNew = 1;
          // Update vertice halfedges          
          v0->halfedge() = h10;
          v1->halfedge() = h1;
          v2->halfedge() = h12;
          v3->halfedge() = h5;
          v->halfedge() = h0;
          
          // Edges
          e0->halfedge() = h0;
          e1->halfedge() = h1;
          e2->halfedge() = h12;
          e3->halfedge() = h14;
          e4->halfedge() = h5;
          e5->halfedge() = h2;
          e6->halfedge() = h10;
          e7->halfedge() = h4;
          
          // Set isNew to true for edges touching the new vertex 
          e0->isNew = 0; // half of original edge
          e6->isNew = 0; // half of original edge
          e5->isNew = 1; // bisecting edge half
          e7->isNew = 1; // bisecting edge half

          // Faces
          f1->halfedge() = h0;
          f2->halfedge() = h3;
          f3->halfedge() = h10;
          f4->halfedge() = h13;
          
          // Return the new vertex
          return v;

      }

    return VertexIter();
  }


    void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
            
        // Set all isNew fields for edges in mesh to false and set newPosition for new vertices
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd();e++) {
            // Set isNew to false for all old edges
            e->isNew = 0;
            
            
            // Inner Halfedges
            HalfedgeIter h0 = e->halfedge();
            HalfedgeIter h1 = h0->next();
            HalfedgeIter h2 = h1->next();
            HalfedgeIter h3 = h0->twin();
            HalfedgeIter h4 = h3->next();
            HalfedgeIter h5 = h4->next();

            // Vertices
            VertexIter v0 = h0->vertex();  // first vertex at endpoint of previous edge
            VertexIter v1 = h3->vertex();  // second vertex at opposite endpoint of previous edge
            VertexIter v2 = h2->vertex();  // first vertex at endpoint of flipped edge
            VertexIter v3 = h5->vertex();  // second vertex at opposite endpoint of flipped edge

            // Precompute and store position of new vertices in newPosition field of old edges
            e->newPosition = 3.0 / 8.0 * (v0->position + v1->position) + 1.0 / 8.0 * (v2->position + v3->position);
            
        }
        
        // Set all vertex isNew members in old mesh to false and set newPosition for old vertices
        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) { 
            // set isNew to false for all old vertices
            v->isNew = 0;
            
            
            // Set newPosition of old vertices as weighted combination of surrounding old vertices
            HalfedgeIter h = v->halfedge();
            Vector3D original_neighbor_position_sum(0,0,0);
            do {
                original_neighbor_position_sum += h->twin()->vertex()->position;
                h = h->twin()->next();

            } while (h != v->halfedge());
            float n = (float) v->degree();
            float u = (n == 3.0) ? (3.0 / 16.0) : (3.0 / (8.0 * n));
            
            // Set newPostion field for old vertices
            v->newPosition = (1.0 - n * u) * v->position + u * original_neighbor_position_sum;
            
        }
        
        // Split all old edges 
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            VertexIter v1 = e->halfedge()->vertex();
            VertexIter v2 = e->halfedge()->twin()->vertex();

            // Split only if its an old edge (prevents infinite loop)
            if(!(v1->isNew||v2->isNew)) {
                VertexIter v = mesh.splitEdge(e);
                v->newPosition = e->newPosition; // set new vertex position
            }
        }

        // Flip appropriate new edges (new edges that have exactly one new vertex)
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            VertexIter v1 = e->halfedge()->vertex();
            VertexIter v2 = e->halfedge()->twin()->vertex();

            if(e->isNew && (v1->isNew + v2->isNew == 1)) {
                mesh.flipEdge(e);  // automatically updates newPosition field in new vertices
            }
        }
        
        
        // Set all vertex positions to their new positions after subdivision and revert isNew to false
        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->position = v->newPosition;
            v->isNew = 0;
        }    
    }
}