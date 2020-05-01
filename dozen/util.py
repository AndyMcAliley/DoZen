# DoZen utilities

# import matplotlib
# matplotlib.use("TkAgg")
# import matplotlib.pyplot as plt
import math
import tkinter as tk
from tkinter import filedialog


def getdir(initialdir = './',title='Select directory',mustexist=True):
    # get root directory from user
    tkroot = tk.Tk()
    tkroot.withdraw()
    # Bring file dialog to front
    # Make it almost invisible - no decorations, 0 size, top left corner.
    tkroot.overrideredirect(True)
    tkroot.geometry('0x0+0+0')
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    tkroot.deiconify()
    tkroot.lift()
    tkroot.focus_force()
    # Call the dialog!
    # rootdir = filedialog.askdirectory()
    rootdir = filedialog.askdirectory(
            initialdir=initialdir,
            title=title,
            mustexist=mustexist)
    # ensure that window closes immediately after mouse click
    tkroot.update()
    # tkroot.quit()
    tkroot.destroy()
    return rootdir


def getfile(initialdir = './',title='Select file'):
    # get root directory from user
    tkroot = tk.Tk()
    tkroot.withdraw()
    # Bring file dialog to front
    # Make it almost invisible - no decorations, 0 size, top left corner.
    tkroot.overrideredirect(True)
    tkroot.geometry('0x0+0+0')
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    tkroot.deiconify()
    tkroot.lift()
    tkroot.focus_force()
    # Call the dialog!
    filename = filedialog.askopenfilename(
            initialdir=initialdir,
            title=title)
    # ensure that window closes immediately after mouse click
    tkroot.update()
    # tkroot.quit()
    tkroot.destroy()
    return filename


def savefile(initialdir = './',title='Select file'):
    # get root directory from user
    tkroot = tk.Tk()
    tkroot.withdraw()

    # tkroot.call('wm', 'attributes', '.', '-topmost', True)

    # Bring file dialog to front
    # Make it almost invisible - no decorations, 0 size, top left corner.
    tkroot.overrideredirect(True)
    tkroot.geometry('0x0+0+0')
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    tkroot.deiconify()
    tkroot.lift()
    tkroot.focus_force()

    # Call the dialog!
    filename = filedialog.asksaveasfilename(
        title=title,
        initialdir=initialdir)
    # ensure that window closes immediately after mouse click
    tkroot.update()
    # tkroot.quit()
    tkroot.destroy()
    return filename


def radian_conversion_factor(units='mrad'):
    '''
    Return scaling factor to convert angle units to radians
    '''
    if (units=='mrad' or units=='milliradians'):
        scale = 0.001
    elif (units == 'rad' or units=='radians'):
        scale = 1.
    elif (units == 'deg' or units=='degrees'):
        scale = math.pi/180
    else:
        raise ValueError('invalid units')
    return scale


def dropdown(names):
    tkroot = tk.Tk()
    # tkroot.withdraw()

    # print('1')
    # names = campaign_settings.name.tolist()
    name = names[0]

    # print('2')
    variable = tk.StringVar(tkroot)
    variable.set(name) # default value

    # print('3')
    # w = apply(tk.OptionMenu, (tkroot, variable) + tuple(names))
    w = tk.OptionMenu(tkroot, variable, *names)
    w.pack()

    # print('4')
    def ok():
        name = variable.get()
        print ("value is " + name)
        tkroot.quit()
        # return name

    # print('5')
    button = tk.Button(tkroot, text="OK", command=ok)
    button.pack()

    # print('6')
    tkroot.call('wm', 'attributes', '.', '-topmost', True)
    t=tk.mainloop()

    # print('7')
    # print(name)
    # ensure that window closes immediately after mouse click
    tkroot.update()
    tkroot.destroy()
    return name


