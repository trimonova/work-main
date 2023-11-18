/*
	This file was generated by XVT-Design 4.5, a product of:

		XVT Software Inc.
		4900 Pearl East Circle
		Boulder, CO USA 80301
		303-443-4223, fax 303-443-0969

	Generated on Tue Sep 02 10:35:00 1997
*/
/* TAG BEGIN SPCL:Pre_Header */
/* TAG END SPCL:Pre_Header */

#include "xvt.h"
#include "xvtcm.h"
#include "nodInc.h"

/* 
	Information about the window
*/
#define WIN_RES_ID BACKGROUND_WINDOW
#define WIN_FLAGS 0x810L
#define WIN_CLASS ""
#define WIN_BORDER W_DOC
/* TAG BEGIN SPCL:Obj_Decl */
/* TAG END SPCL:Obj_Decl */

/*
	Handler for window BACKGROUND_WINDOW ("Background Color")
*/
long XVT_CALLCONV1
#if XVT_CC_PROTO
BACKGROUND_WINDOW_eh XVT_CALLCONV2 (WINDOW xdWindow, EVENT *xdEvent)
#else
BACKGROUND_WINDOW_eh XVT_CALLCONV2 (xdWindow, xdEvent)
WINDOW xdWindow;
EVENT *xdEvent;
#endif
{
	short xdControlId = xdEvent->v.ctl.id;
	/* TAG BEGIN SPCL:Var_Decl */
	/* TAG END SPCL:Var_Decl */

	switch (xdEvent->type) {
	case E_CREATE:
		/*
			Window has been created; first event sent to newly-created
			window.
		*/
		{
		/* TAG BEGIN EVNT:Create */
		/* TAG END EVNT:Create */
		}
		break;
	case E_DESTROY:
		/*
			Window has been closed; last event sent to window.
		*/
		xdRemoveHelpAssoc( xdWindow );
		{
		/* TAG BEGIN EVNT:Destroy */
		/* TAG END EVNT:Destroy */
		}
		break;
	case E_FOCUS:
		{
		/* TAG BEGIN EVNT:Focus */
		/*
			Window has lost or gained focus.
		*/
		if (xdEvent->v.active) {
			/*
				Window has gained focus
			*/
		} else {
			/*
				Window has lost focus
			*/
		}
		/* TAG END EVNT:Focus */
		}
		break;
	case E_SIZE:
		/*
			Size of window has been set or changed; sent when window is
			created or subsequently resized by user or via xvt_vobj_move.
		*/
		{
		/* TAG BEGIN EVNT:Size */
		/* TAG END EVNT:Size */
		}
		break;
	case E_UPDATE:
		/*
			Window requires updating.
		*/
		{
		/* TAG BEGIN EVNT:Update */
		xvt_dwin_clear(xdWindow, (COLOR)xvt_vobj_get_attr(xdWindow, ATTR_BACK_COLOR));
		/* TAG END EVNT:Update */
		}
		break;
	case E_CLOSE:
		/*
			Request to close window; user operated "close" menu item on
			window system menu, or operated "close" control on window
			frame. Not sent if Close on File menu is issued. Window not
			closed unless xvt_vobj_destroy is called.
		*/
		{
		/* TAG BEGIN EVNT:Close */
		xvt_vobj_destroy(xdWindow);
		/* TAG END EVNT:Close */
		}
		break;
	case E_CHAR:
		/*
			Character typed.
		*/
		{
		/* TAG BEGIN EVNT:Char */
		/* TAG END EVNT:Char */
		}
		break;
	case E_MOUSE_UP:
		/*
			Mouse was released
		*/
		{
		/* TAG BEGIN EVNT:Mouse_Up */
		/* TAG END EVNT:Mouse_Up */
		}
		break;
	case E_MOUSE_DOWN:
		/*
			Mouse was pressed
		*/
		{
		/* TAG BEGIN EVNT:Mouse_Down */
		/* TAG END EVNT:Mouse_Down */
		}
		break;
	case E_MOUSE_DBL:
		/*
			Mouse was double clicked
		*/
		{
		/* TAG BEGIN EVNT:Mouse_Dbl */
		/* TAG END EVNT:Mouse_Dbl */
		}
		break;
	case E_MOUSE_MOVE:
		/*
			Mouse was moved
		*/
		{
		/* TAG BEGIN EVNT:Mouse_Move */
		/* TAG END EVNT:Mouse_Move */
		}
		break;
	case E_HSCROLL:
		{
		/* TAG BEGIN EVNT:Hscroll */
		/*
			Horizontal scrollbar on frame was operated
		*/
		switch (xdEvent->v.scroll.what) {
		case SC_LINE_UP:
			break;
		case SC_LINE_DOWN:
			break;
		case SC_PAGE_UP:
			break;
		case SC_PAGE_DOWN:
			break;
		case SC_THUMB:
			break;
		case SC_THUMBTRACK:
			break;
		default:
			break;
		}
		/* TAG END EVNT:Hscroll */
		}
		break;
	case E_VSCROLL:
		{
		/* TAG BEGIN EVNT:Vscroll */
		/*
			Vertical scrollbar on frame was operated
		*/
		switch (xdEvent->v.scroll.what) {
		case SC_LINE_UP:
			break;
		case SC_LINE_DOWN:
			break;
		case SC_PAGE_UP:
			break;
		case SC_PAGE_DOWN:
			break;
		case SC_THUMB:
			break;
		case SC_THUMBTRACK:
			break;
		default:
			break;
		}
		/* TAG END EVNT:Vscroll */
		}
		break;
	case E_COMMAND:
		/*
			User issued command on window menu bar (menu bar at top of
			screen for Mac/CH).
		*/
		{
		/* TAG BEGIN EVNT:Command */
		/*
			No menubar was associated with this window
		*/
		/* TAG END EVNT:Command */
		}
		break;
	case E_CONTROL:
		/*
			User operated control in window.
		*/
		{
		/* TAG BEGIN SPCL:Control_Decl */
		/* TAG END SPCL:Control_Decl */

		switch(xdControlId) {
		case BACKGROUND_BUTTON: /* "Background Color..." */
			{
			/* TAG BEGIN BACKGROUND_BUTTON EVNT:Control */
			WINDOW win;
			if (!(win = createCenteredWindow(DEFINE_COLOUR_WINDOW, TASK_WIN, EM_ALL,
			                                                      DEFINE_COLOUR_WINDOW_eh, 0L)))
				xvt_dm_post_error("Can't open window");
			/* TAG END BACKGROUND_BUTTON EVNT:Control */
			}
			break;
		default:
			break;
		}
		}
		break;
	case E_FONT:
		/*
			User issued font command on window menu bar (menu bar at top of
			screen for Mac/CH).
		*/
		{
		/* TAG BEGIN EVNT:Font */
		/* TAG END EVNT:Font */
		}
		break;
	case E_TIMER:
		/*
			Timer associated with window went off.
		*/
		{
		/* TAG BEGIN EVNT:Timer */
		/* TAG END EVNT:Timer */
		}
		break;
	case E_USER:
		/*
			Application initiated.
		*/
		{
		/* TAG BEGIN EVNT:User */
		switch (xdEvent->v.user.id) {
		case -1:
		default:
			break;
		}
		/* TAG END EVNT:User */
		}
		break;
	default:
		break;
	}
	/* TAG BEGIN SPCL:Bottom */
#ifdef XVT_R40_TXEDIT_BEHAVIOR
	xvt_tx_process_event(xdWindow, xdEvent);
#endif
	/* TAG END SPCL:Bottom */
	return 0L;
}
