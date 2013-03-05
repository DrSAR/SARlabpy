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

class Scan ( HasTraits ):
    name  = Str( '<unknown>' )
    title = Str
    phone = Regex( regex = r'\d\d\d-\d\d\d\d' )

    def default_title ( self ):
        self.title = 'Senior Engineer'

class Scanlist ( HasTraits ):
    name        = Str( '<unknown>' )
    scanlist   = List( Scan )

class Owner ( HasTraits ):
    name    = Str( '<unknown>' )
    scanlist = Instance( Scanlist )

# INSTANCES

jason = Scan(
     name  = 'Jason',
     title = 'Engineer',
     phone = '536-1057' )

mike = Scan(
     name  = 'Mike',
     title = 'Sr. Marketing Analyst',
     phone = '536-1057' )

dave = Scan(
     name  = 'Dave',
     title = 'Sr. Engineer',
     phone = '536-1057' )

susan = Scan(
     name  = 'Susan',
     title = 'Engineer',
     phone = '536-1057' )

betty = Scan(
     name  = 'Betty',
     title = 'Marketing Analyst' )

owner = Owner(
    name    = 'wile',
    scanlist = Scanlist(
        name = 'Acme Labs, Inc.',
        scanlist = [ dave, susan, mike, betty, jason ]
    )
)

# View for objects that aren't edited
no_view = View()

# Actions used by tree editor context menu

def_title_action = Action(name='Default title',
                          action = 'object.default')

# View used by tree editor
scanlist_view = View(
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
        TreeNode( node_for  = [ Scanlist ],
                  auto_open = False,
                  children  = 'scanlist',
                  label     = '=Scans',
                  view      = no_view,
                  add       = [ Scan ] ),
        TreeNode( node_for  = [ Scan ],
                  auto_open = False,
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
                  view = scanlist_view )
    ]
)

# The main view
view = View(
           Group(
               Item(
                    name = 'scanlist',
                    id = 'scanlist',
                    editor = tree_editor,
                    resizable = True ),
                orientation = 'vertical',
                show_labels = True,
                show_left = True, ),
            title = 'View of Scans',
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