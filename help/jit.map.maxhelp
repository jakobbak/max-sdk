{
    "patcher": {
        "fileversion": 1,
        "appversion": {
            "major": 7,
            "minor": 0,
            "revision": 0,
            "architecture": "x86",
            "modernui": 1
        },
        "rect": [100, 100, 535, 482],
        "bglocked": 0,
        "openinpresentation": 0,
        "default_fontsize": 12,
        "default_fontface": 0,
        "default_fontname": "Arial",
        "gridonopen": 0,
        "gridsize": [15, 15],
        "gridsnaponopen": 0,
        "statusbarvisible": 2,
        "toolbarvisible": 1,
        "lefttoolbarpinned": 0,
        "toptoolbarpinned": 0,
        "righttoolbarpinned": 0,
        "bottomtoolbarpinned": 0,
        "toolbars_unpinned_last_save": 0,
        "tallnewobj": 0,
        "boxanimatetime": 200,
        "enablehscroll": 1,
        "enablevscroll": 1,
        "devicewidth": 0,
        "description": "",
        "digest": "",
        "tags": "",
        "style": "",
        "showrootpatcherontab": 0,
        "showontab": 0,
        "boxes": [
            {
                "box": {
                    "fontname": "Arial",
                    "fontsize": 12,
                    "id": "obj-5",
                    "maxclass": "newobj",
                    "numinlets": 1,
                    "numoutlets": 1,
                    "outlettype": [""],
                    "patching_rect": [125, 139, 133, 20],
                    "saved_object_attributes": {
                        "filename": "helpstarter.js",
                        "parameter_enable": 0
                    },
                    "style": "",
                    "text": "js helpstarter.js jit.map"
                }
            },
            {
                "box": {
                    "fontname": "Arial",
                    "fontsize": 12,
                    "id": "obj-6",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 7,
                            "minor": 0,
                            "revision": 0,
                            "architecture": "x86",
                            "modernui": 1
                        },
                        "rect": [100, 126, 535, 456],
                        "bglocked": 0,
                        "openinpresentation": 0,
                        "default_fontsize": 13,
                        "default_fontface": 0,
                        "default_fontname": "Arial",
                        "gridonopen": 0,
                        "gridsize": [15, 15],
                        "gridsnaponopen": 0,
                        "statusbarvisible": 2,
                        "toolbarvisible": 1,
                        "lefttoolbarpinned": 0,
                        "toptoolbarpinned": 0,
                        "righttoolbarpinned": 0,
                        "bottomtoolbarpinned": 0,
                        "toolbars_unpinned_last_save": 0,
                        "tallnewobj": 0,
                        "boxanimatetime": 200,
                        "enablehscroll": 1,
                        "enablevscroll": 1,
                        "devicewidth": 0,
                        "description": "",
                        "digest": "",
                        "tags": "",
                        "style": "",
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-3",
                                    "linecount": 2,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [214, 133, 230, 36],
                                    "style": "helpfile_label",
                                    "text": "map values are always float numbers regardless of matrix input type"
                                }
                            },
                            {
                                "box": {
                                    "colwidth": 40,
                                    "fontface": 0,
                                    "fontname": "Arial",
                                    "fontsize": 13,
                                    "hscroll": 0,
                                    "id": "obj-10",
                                    "maxclass": "jit.cellblock",
                                    "numinlets": 2,
                                    "numoutlets": 4,
                                    "outlettype": ["list", "", "", ""],
                                    "patching_rect": [185, 285, 321, 144],
                                    "vscroll": 0
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-8",
                                    "maxclass": "jit.pwindow",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": ["", ""],
                                    "patching_rect": [28, 285, 148, 144],
                                    "presentation_rect": [30, 30, 128, 128]
                                }
                            },
                            {
                                "box": {
                                    "attr": "clip",
                                    "displaymode": 8,
                                    "fontface": 0,
                                    "fontname": "Arial",
                                    "fontsize": 13,
                                    "id": "obj-13",
                                    "lock": 1,
                                    "maxclass": "attrui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [""],
                                    "patching_rect": [247, 213, 74, 23],
                                    "style": "",
                                    "text_width": 47
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-11",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": ["bang"],
                                    "patching_rect": [14, 162, 24, 24],
                                    "style": ""
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-9",
                                    "maxclass": "toggle",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": ["int"],
                                    "parameter_enable": 0,
                                    "patching_rect": [41, 134, 24, 24],
                                    "style": ""
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 13,
                                    "id": "obj-7",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": ["bang"],
                                    "patching_rect": [41, 162, 71, 23],
                                    "style": "",
                                    "text": "qmetro 33"
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 13,
                                    "id": "obj-6",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": ["jit_matrix", ""],
                                    "patching_rect": [41, 199, 148, 23],
                                    "style": "",
                                    "text": "jit.noise 1 float32 10 10"
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 13,
                                    "id": "obj-5",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": ["jit_matrix", ""],
                                    "patching_rect": [185, 248, 50, 23],
                                    "style": "",
                                    "text": "jit.map"
                                }
                            },
                            {
                                "box": {
                                    "border": 0,
                                    "filename": "helpdetails.js",
                                    "id": "obj-2",
                                    "ignoreclick": 1,
                                    "jsarguments": ["jit.map", 85],
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [""],
                                    "parameter_enable": 0,
                                    "patching_rect": [6, 5, 500, 118]
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [1, 0.788235, 0.470588, 1],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-93",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": ["", "", "int"],
                                    "parameter_enable": 0,
                                    "patching_rect": [66.5, 134, 20, 20],
                                    "rounded": 60,
                                    "style": "",
                                    "text": "1",
                                    "textcolor": [0.34902, 0.34902, 0.34902, 1]
                                }
                            },
                            {
                                "box": {
                                    "attr": "map",
                                    "displaymode": 4,
                                    "fontface": 0,
                                    "fontname": "Arial",
                                    "fontsize": 13,
                                    "id": "obj-12",
                                    "maxclass": "attrui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [""],
                                    "patching_rect": [214, 179, 230, 23],
                                    "style": "",
                                    "text_width": 61
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": ["obj-6", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-11", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-5", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-12", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-5", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-13", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-10", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-5", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-8", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-5", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-5", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-6", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-6", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-7", 0]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": ["obj-7", 0],
                                    "disabled": 0,
                                    "hidden": 0,
                                    "source": ["obj-9", 0]
                                }
                            }
                        ]
                    },
                    "patching_rect": [10, 85, 50, 20],
                    "saved_object_attributes": {
                        "description": "",
                        "digest": "",
                        "fontsize": 13,
                        "globalpatchername": "",
                        "style": "",
                        "tags": ""
                    },
                    "style": "",
                    "text": "p basic",
                    "varname": "basic_tab"
                }
            },
            {
                "box": {
                    "border": 0,
                    "filename": "helpname.js",
                    "id": "obj-8",
                    "ignoreclick": 1,
                    "jsarguments": ["jit.map"],
                    "maxclass": "jsui",
                    "numinlets": 1,
                    "numoutlets": 1,
                    "outlettype": [""],
                    "parameter_enable": 0,
                    "patching_rect": [10, 10, 305, 55]
                }
            },
            {
                "box": {
                    "fontname": "Arial",
                    "fontsize": 12,
                    "id": "obj-1",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 7,
                            "minor": 0,
                            "revision": 0,
                            "architecture": "x86",
                            "modernui": 1
                        },
                        "rect": [0, 26, 535, 456],
                        "bglocked": 0,
                        "openinpresentation": 0,
                        "default_fontsize": 13,
                        "default_fontface": 0,
                        "default_fontname": "Arial",
                        "gridonopen": 0,
                        "gridsize": [15, 15],
                        "gridsnaponopen": 0,
                        "statusbarvisible": 2,
                        "toolbarvisible": 1,
                        "lefttoolbarpinned": 0,
                        "toptoolbarpinned": 0,
                        "righttoolbarpinned": 0,
                        "bottomtoolbarpinned": 0,
                        "toolbars_unpinned_last_save": 0,
                        "tallnewobj": 0,
                        "boxanimatetime": 200,
                        "enablehscroll": 1,
                        "enablevscroll": 1,
                        "devicewidth": 0,
                        "description": "",
                        "digest": "",
                        "tags": "",
                        "style": "",
                        "showontab": 1,
                        "boxes": [],
                        "lines": []
                    },
                    "patching_rect": [205, 205, 50, 20],
                    "saved_object_attributes": {
                        "description": "",
                        "digest": "",
                        "fontsize": 13,
                        "globalpatchername": "",
                        "style": "",
                        "tags": ""
                    },
                    "style": "",
                    "text": "p ?",
                    "varname": "q_tab"
                }
            }
        ],
        "lines": [],
        "dependency_cache": [
            {
                "name": "helpname.js",
                "bootpath": "/Applications/Max 7/Max.app/Contents/Resources/C74/help/resources",
                "patcherrelativepath": "../../../../../../Applications/Max 7/Max.app/Contents/Resources/C74/help/resources",
                "type": "TEXT",
                "implicit": 1
            },
            {
                "name": "helpdetails.js",
                "bootpath": "/Applications/Max 7/Max.app/Contents/Resources/C74/help/resources",
                "patcherrelativepath": "../../../../../../Applications/Max 7/Max.app/Contents/Resources/C74/help/resources",
                "type": "TEXT",
                "implicit": 1
            },
            {
                "name": "helpstarter.js",
                "bootpath": "/Applications/Max 7/Max.app/Contents/Resources/C74/help/resources",
                "patcherrelativepath": "../../../../../../Applications/Max 7/Max.app/Contents/Resources/C74/help/resources",
                "type": "TEXT",
                "implicit": 1
            }
        ],
        "embedsnapshot": 0
    }
}