using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using Plankton;
using PlanktonGh;

/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(Mesh mesh, List<Point3d> ct, ref object verticeAdjFacesId, ref object verticeIdAtAdjFaces, ref object planktonMesh, ref object averagePlane, ref object faceVIndices, ref object adjFaceIndices, ref object vertices, ref object thicknessDomains, ref object periPlanes, ref object cenPlanes)
  {
    var dual = mesh.ToPlanktonMesh().Dual();
    dual.ReplaceVertices(ct);

    DataTree<int> vf_adj = new DataTree<int>();
    DataTree<int> vf_adj_id = new DataTree<int>();

    var vs = dual.Vertices;
    var hes = dual.Halfedges;

    List<Plane> pls = new List<Plane>();
    List<Interval> domains = new List<Interval>();
    DataTree<int> adjf_indices = new DataTree<int>();
    DataTree<int> pts_indices = new DataTree<int>();


    Point3d[] vts = new Point3d[dual.Vertices.Count];
    for(int i = 0;i < dual.Vertices.Count;i++)
    {
      vts[i] = vs[i].ToPoint3d();
      List<int> f_id = new List<int>();
      List<int> v_at_f_id = new List<int>();
      foreach(int he in hes.GetVertexCirculator(vs[i].OutgoingHalfedge))
      {
        int face = dual.Halfedges[he].AdjacentFace;
        if(face >= 0)
        {
          f_id.Add(face);
          int id = 0;
          foreach(int she in dual.Faces.GetHalfedgesCirculator(face))
          {
            if(hes[she].StartVertex == i)break;
            id++;
          }
          v_at_f_id.Add(id);
        }
      }
      vf_adj_id.AddRange(v_at_f_id, new GH_Path(i));
      vf_adj.AddRange(f_id, new GH_Path(i));
    }
    verticeAdjFacesId = vf_adj_id;
    verticeIdAtAdjFaces = vf_adj;

    for(int i = 0;i < dual.Faces.Count;i++)
    {
      Plane pl;
      Point3d pA = new Point3d();
      Interval dom = new Interval(double.MaxValue, double.MinValue);
      List<Point3d> pts = new List<Point3d>();
      int count = 0;
      foreach(int he_index in dual.Faces.GetHalfedgesCirculator(i))
      {
        Point3d pt = dual.Vertices[dual.Halfedges[he_index].StartVertex].ToPoint3d();
        pts.Add(pt); pA += pt;
        pts_indices.Add(dual.Halfedges[he_index].StartVertex, new GH_Path(i));
        adjf_indices.Add(dual.Halfedges[dual.Halfedges.GetPairHalfedge(he_index)].AdjacentFace, new GH_Path(i));
        count++;
      }
      pA /= count;
      Plane.FitPlaneToPoints(pts, out pl);
      pl.Origin = pA;

      Vector3d vA = Vector3d.Zero;
      for(int k = 0;k < pts.Count;k++)
      {
        Point3d rpt = Point3d.Origin;
        pl.RemapToPlaneSpace(pts[k], out rpt);
        dom.T1 = dom.T1 > rpt.Z ? dom.T1 : rpt.Z;
        dom.T0 = dom.T0 < rpt.Z ? dom.T0 : rpt.Z;
        dom.T1 = dom.T1 > 0.05 ? dom.T1 : 0.05;
        dom.T0 = dom.T0 < -0.05 ? dom.T0 : -0.05;
        vA += Vector3d.CrossProduct(pts[k] - pA, pts[k < (pts.Count - 1) ? (k + 1) : 0] - pA);
      }
      vA.Unitize();
      if (vA * pl.ZAxis < 0)
      {
        pl.Flip();dom.T1 = 0 - dom.T0;dom.T0 = 0 - dom.T1;
      }
      pls.Add(pl);
      domains.Add(dom);
    }

    planktonMesh = dual;
    averagePlane = pls;
    faceVIndices = pts_indices;
    adjFaceIndices = adjf_indices;
    vertices = vts;
    thicknessDomains = domains;

    //construct solid
    DataTree<Plane> peri_planes = new DataTree<Plane>();
    DataTree<Plane> cen_planes = new DataTree<Plane>();
    for(int f_index = 0;f_index < dual.Faces.Count;f_index++)
    {
      int adjf_index,ptStart_index,ptEnd_index;
      int num = adjf_indices.Branch(pts_indices.Path(f_index)).Count;
      List<Plane> peri_pls = new  List<Plane>();
      for(int j = 0;j < num;j++)
      {
        ptStart_index = pts_indices[pts_indices.Path(f_index), j];
        ptEnd_index = pts_indices[pts_indices.Path(f_index), j == (num - 1) ? 0 : j + 1];
        Point3d s = vts[ptStart_index];
        Point3d e = vts[ptEnd_index];

        adjf_index = adjf_indices[adjf_indices.Path(f_index), j];
        Vector3d vecY;
        if(adjf_index != -1)
          vecY = pls[f_index].ZAxis + pls[adjf_index].ZAxis;
        else
          vecY = pls[f_index].ZAxis;
        peri_pls.Add(new Plane((s + e) / 2, s - e, vecY));
      }
      peri_planes.AddRange(peri_pls, new GH_Path(f_index));


      Plane[] cen_pls = new Plane[2]{pls[f_index],pls[f_index]};
      cen_pls[1].Translate(cen_pls[1].ZAxis * domains[f_index].T1);
      cen_pls[0].Translate(cen_pls[0].ZAxis * domains[f_index].T0);
      cen_pls[1].Flip();
      cen_planes.AddRange(cen_pls, new GH_Path(f_index));
    }
    periPlanes = peri_planes;
    cenPlanes = cen_planes;
  }

  // <Custom additional code> 

  // </Custom additional code> 
}
