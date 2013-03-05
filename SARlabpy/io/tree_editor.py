# -*- coding: utf-8 -*-
"""
`Example 18: Code for example tree editor <http://code.enthought.com/projects/traits/docs/html/TUIUG/factories_advanced_extra.html#treeeditor>`_
"""

# tree_editor.py -- Example of a tree editor

from enthought.traits.api \
    import HasTraits, Str, Regex, List, Instance
from enthought.traits.ui.api \
    import TreeEditor, TreeNode, View, Item, VSplit, \
           HGroup, Handler, Group
from enthought.traits.ui.menu \
    import Menu, Action, Separator
from enthought.traits.ui.wx.tree_editor \
    import NewAction, CopyAction, CutAction, \
           PasteAction, DeleteAction, RenameAction

# DATA CLASSES

class Employee ( HasTraits ):
    name  = Str( '<unknown>' )
    title = Str
    phone = Regex( regex = r'\d\d\d-\d\d\d\d' )

    def default_title ( self ):
        self.title = 'Senior Engineer'

class Company ( HasTraits ):
    name        = Str( '<unknown>' )
    employees   = List( Employee )

class Owner ( HasTraits ):
    name    = Str( '<unknown>' )
    company = Instance( Company )

# INSTANCES

jason = Employee(
     name  = 'Jason',
     title = 'Engineer',
     phone = '536-1057' )

mike = Employee(
     name  = 'Mike',
     title = 'Sr. Marketing Analyst',
     phone = '536-1057' )

dave = Employee(
     name  = 'Dave',
     title = 'Sr. Engineer',
     phone = '536-1057' )

susan = Employee(
     name  = 'Susan',
     title = 'Engineer',
     phone = '536-1057' )

betty = Employee(
     name  = 'Betty',
     title = 'Marketing Analyst' )

owner = Owner(
    name    = 'wile',
    company = Company(
        name = 'Acme Labs, Inc.',
        employees = [ dave, susan, mike, betty, jason ]
    )
)

# View for objects that aren't edited
no_view = View()

# Actions used by tree editor context menu

def_title_action = Action(name='Default title',
                          action = 'object.default')

# View used by tree editor
employee_view = View(
    VSplit(
        HGroup( '3', 'name' ),
        HGroup( '9', 'title' ),
        HGroup( 'phone' ),
        id = 'vsplit' ),
    id = 'enthought.traits.doc.example.treeeditor',
    dock = 'vertical' )

# Tree editor
tree_editor = TreeEditor(
    nodes = [
        TreeNode( node_for  = [ Company ],
                  auto_open = True,
                  children  = 'employees',
                  label     = '=Employees',
                  view      = no_view,
                  add       = [ Employee ] ),
        TreeNode( node_for  = [ Employee ],
                  auto_open = True,
                  label     = 'name',
                  menu=Menu( NewAction,
                             Separator(),
                             def_title_action,
                             Separator(),
                             CopyAction,
                             CutAction,
                             PasteAction,
                             Separator(),
                             DeleteAction,
                             Separator(),
                             RenameAction ),
                  view = employee_view )
    ]
)

# The main view
view = View(
           Group(
               Item(
                    name = 'company',
                    id = 'company',
                    editor = tree_editor,
                    resizable = True ),
                orientation = 'vertical',
                show_labels = True,
                show_left = True, ),
            title = 'Company Structure',
            id = \
             'enthought.traits.ui.tests.tree_editor_test',
            dock = 'horizontal',
            drop_class = HasTraits,
            buttons = [ 'Undo', 'OK', 'Cancel' ],
            resizable = True,
            width = .3,
            height = .3 )

if __name__ == '__main__':
    owner.configure_traits( view = view )