# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 00:07:53 2024

@author: valerie
"""

from PySide6.QtWidgets import QMainWindow, QApplication
from PySide6.QtCore import QFile, QIODevice
import pyqtgraph as pg
import sys
from PySide6.QtUiTools import QUiLoader
from PySide6 import QtGui
from ui_mainwindow import Ui_MainWindow
from pyqtgraph import PlotWidget


class MainWindow(QMainWindow):
    
    def __init__(self):
        super(MainWindow, self).__init__()
        
        ui_file2 =True
        if ui_file2 == True :
            print("ui file")
            loader = QUiLoader()
            loader.registerCustomWidget(PlotWidget)
            ui_file_name = "img.ui"
            ui_file = QFile(ui_file_name)
            
            
            if not ui_file.open(QIODevice.ReadOnly):
                print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
                sys.exit(-1)
            
            self.ui = loader.load(ui_file)
            ui_file.close()
        else : 
            self.ui = Ui_MainWindow()
            self.ui.setupUi(self)

        #sys.stdout = Log(self.ui.log_edit)
        #self.graphWidget = pg.PlotWidget()
        #self.setCentralWidget(self.graphWidget)

        hour = [1,2,3,4,5,6,7,8,9,10]
        temperature = [30,32,34,32,33,31,29,32,35,45]

        self.ui.graphWidget.setBackground('w')
        self.ui.graphWidget.plot(hour, temperature)
        self.ui.print_btn.clicked.connect (self.print_clicked)
        #print("Message console on textEdit")
        
        app.aboutToQuit.connect(self.closeEvent)
    
    def print_clicked(self) :
        print('message print')
    
    def show(self) :
        self.ui.show()
        
    def closeEvent (self):
        print('fire ')
        
        
"""
class Log(object):
    def __init__(self, edit):
        self.out = sys.stdout
        self.textEdit = edit

    def write(self, message):
        self.out.write(message)
        #self.textEdit.append(message)
        self.textEdit.moveCursor(QtGui.QTextCursor.End) 
        self.textEdit.insertPlainText( message )

    def flush(self):
        self.out.flush()
"""

if __name__ == "__main__":
    
    ui_file= True
    
    ## doit etre initialiser ici
    if ui_file==True :
        loader = QUiLoader()
    
    
    # pour eviter de devoir tuer app qt sous spyder
    app = QApplication.instance() 
    if not app:
        app = QApplication(sys.argv)
    else:
        app = QApplication.instance()
    
    my_wnd_class=MainWindow() 
    
    my_wnd_class.show()

    sys.stdout.flush()
    sys.exit(app.exec())