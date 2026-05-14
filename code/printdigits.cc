#include "DaqEventsManager.h"
#include "ChipSADC.h"
#include <iostream>
#include <set>
using namespace CS;
using namespace std;

int main(int argc, char *argv[]) {
  DaqEventsManager manager;
  manager.SetMapsDir("../maps");
  manager.AddDataSource(argv[1]);

  int nevt = 0;
  while (manager.ReadEvent() && nevt < 10) {
    if (!manager.DecodeEvent()) continue;
    nevt++;
    cout << "=== Event " << nevt << " ===" << endl;
    set<string> seen;
    for (auto it = manager.GetEventDigits().begin();
              it != manager.GetEventDigits().end(); ++it) {
      const ChipSADC::Digit* s = dynamic_cast<const ChipSADC::Digit*>(it->second);
      if (s) {
        string n = s->GetDetID().GetName();
        if (!seen.count(n)) {
          cout << "  SADC: " << n << endl;
          seen.insert(n);
        }
      }
    }
  }
  return 0;
}
