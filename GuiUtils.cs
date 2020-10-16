using System;
using System.Text.RegularExpressions;
using UnityEngine;

namespace HopperGuidance
{
    static public class GuiUtils
    {
        public static bool OverPartWindow(Part part, Vector2d mouse)
        {
          bool symmetry = false;
          UIPartActionWindow window = UIPartActionController.Instance.GetItem(part, symmetry);
          return window.Hover;
        }

        public static bool GetMouseHit(CelestialBody body, out RaycastHit hit, Part part)
        {
          if ((part!=null) && OverPartWindow(part, Input.mousePosition))
          {
            hit = new RaycastHit();
            return false;
          }

          // Cast a ray from screen point
          Ray ray = FlightCamera.fetch.mainCamera.ScreenPointToRay(Input.mousePosition);
          bool isHit = Physics.Raycast(ray.origin, ray.direction, out hit, Mathf.Infinity, 1<<15);
          return isHit;
        }
    }
}
