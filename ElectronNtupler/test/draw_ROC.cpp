#include "ROC_drawer.cpp"


void draw_ROC()
{
  ROC_Drawer MyROCDrawer;
  MyROCDrawer.name = "ROC_endcap";
  MyROCDrawer.addSelection = " isEB != 1";
  MyROCDrawer.draw_ROC();
  MyROCDrawer.name = "ROC_barrel";
  MyROCDrawer.addSelection = " isEB == 1";
  MyROCDrawer.draw_ROC();
}