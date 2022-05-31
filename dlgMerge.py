def show_exception_and_exit(exc_type, exc_value, tb):
    import traceback
    traceback.print_exception(exc_type, exc_value, tb)
    input("Press key to exit.")
    sys.exit(-1)

import sys

sys.excepthook = show_exception_and_exit

#Check if any input is not a dlg file
if (len([x for x in sys.argv[1:] if ".dlg" not in x]) != 0):
    raise ValueError("An argument is not a .dlg file.")
if (len(sys.argv) == 1):
    raise ValueError("No .dlg files supplied")

dlgText = ""
for i in range(1,len(sys.argv)):
    dlgFile = open(sys.argv[i],"r")
    dlgLines = dlgFile.read();
    if (i == 1):
        mergedName= [x for x in dlgLines.split("\n") if "DPF> move" in x][0].split(" ")[2].split(".")[0]
    dlgFile.close()
    dlgText += dlgLines;
    dlgText += "\n"
path = "./" + mergedName + "_Merged.dlg"
print(path)
dlgWriter = open(path, "w")
dlgWriter.write(dlgText);

input("Press enter to exit")