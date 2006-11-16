#include "toy-framework.h"
#include "point-ops.h"
#include "point-fns.h"

#include "cairo-features.h"
#if CAIRO_HAS_PDF_SURFACE
#include "cairo-pdf.h"
#endif
#if CAIRO_HAS_SVG_SURFACE
#include "cairo-svg.h"
#endif

GtkWindow* window;
static GtkWidget *canvas;
Toy* current_toy;

//Utility functions

double uniform() {
    return double(rand()) / RAND_MAX;
}

void draw_number(cairo_t *cr, Geom::Point pos, int num) {
    std::ostringstream number;
    number << num;
    cairo_move_to(cr, pos);
    PangoLayout* layout = pango_cairo_create_layout (cr);
      pango_layout_set_text(layout, number.str().c_str(), -1);
      PangoFontDescription *font_desc = pango_font_description_new();
        pango_font_description_set_family(font_desc, "Sans");
      pango_layout_set_font_description(layout, font_desc);
    PangoRectangle logical_extent;
    pango_layout_get_pixel_extents(layout, NULL, &logical_extent);
    pango_cairo_show_layout(cr, layout); 
}

//Framework Accessors
void redraw() { gtk_widget_queue_draw(GTK_WIDGET(window)); }

void Toy::draw(cairo_t *cr, std::ostringstream *notify, int width, int height, bool save)
{
    if(should_draw_bounds()) {
        cairo_set_source_rgba (cr, 0., 0., 0, 0.8);
        cairo_set_line_width (cr, 0.5);
        for(int i = 1; i < 4; i+=2) {
            cairo_move_to(cr, 0, i*width/4);
            cairo_line_to(cr, width, i*width/4);
            cairo_move_to(cr, i*width/4, 0);
            cairo_line_to(cr, i*width/4, height);
        }
    }

    cairo_set_source_rgba (cr, 0., 0.5, 0, 1);
    cairo_set_line_width (cr, 1);
    for(int i = 0; i < handles.size(); i++) {
        draw_circ(cr, handles[i]);
        if(should_draw_numbers()) draw_number(cr, handles[i], i);
    }
    
    cairo_set_source_rgba (cr, 0.5, 0, 0, 1);
    if(selected_handle != NULL && mouse_down == true)
        draw_circ(cr, *selected_handle);

    cairo_set_source_rgba (cr, 0.5, 0.25, 0, 1);
    cairo_stroke(cr);

    cairo_set_source_rgba (cr, 0., 0.5, 0, 0.8);
    {
        *notify << std::ends;
        PangoLayout *layout = gtk_widget_create_pango_layout(GTK_WIDGET(window), notify->str().c_str());
        PangoRectangle logical_extent;
        pango_layout_get_pixel_extents(layout, NULL, &logical_extent);
        cairo_move_to(cr, 0, height-logical_extent.height);
        pango_cairo_show_layout(cr, layout);
    }
}

void Toy::mouse_moved(GdkEventMotion* e)
{
    Geom::Point mouse(e->x, e->y);
    
    if(e->state & (GDK_BUTTON1_MASK | GDK_BUTTON3_MASK)) {
        if(selected_handle != NULL) *selected_handle = mouse;
        redraw();
    }

    if(e->state & (GDK_BUTTON2_MASK)) redraw();
    
    old_mouse_point = mouse;
}

void Toy::mouse_pressed(GdkEventButton* e) {
    Geom::Point mouse(e->x, e->y);
    if(e->button == 1) {
        for(int i = 0; i < handles.size(); i++) {
            if(Geom::L2(mouse - handles[i]) < 5) selected_handle = &handles[i];
        }
        redraw();
        mouse_down = true;
    } else if(e->button == 2) {
        redraw();
    }
    old_mouse_point = mouse;
}

void Toy::mouse_released(GdkEventButton* e) {
    selected_handle = NULL;
    if(e->button == 1) mouse_down = false;
    redraw();
}

//Gui Event Callbacks

void make_about() {
    GtkWidget* about_window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(about_window), "About");
    gtk_window_set_policy(GTK_WINDOW(about_window), FALSE, FALSE, TRUE);
    
    GtkWidget* about_text = gtk_text_view_new();
    GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(about_text));
    gtk_text_buffer_set_text(buf, "Toy lib2geom application", -1);
    gtk_container_add(GTK_CONTAINER(about_window), about_text);

    gtk_widget_show_all(about_window);
}

Geom::Point read_point(FILE* f) {
    Geom::Point p;
    for(int i = 0; i < 2; i++)
        assert(fscanf(f, " %lf ", &p[i]));
    return p;
}

void open() {
    if(current_toy != NULL) {
    GtkWidget* d = gtk_file_chooser_dialog_new("Open handle configuration", window, GTK_FILE_CHOOSER_ACTION_OPEN, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT, NULL);
        if(gtk_dialog_run(GTK_DIALOG(d)) == GTK_RESPONSE_ACCEPT) {
            const char* filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));
            FILE* f = fopen(filename, "r");
            current_toy->handles.clear();
            while(!feof(f))
                current_toy->handles.push_back(read_point(f));
            fclose(f);
        }
        gtk_widget_destroy(d);
    }
}

void save() {
    if(current_toy != NULL) {
        GtkWidget* d = gtk_file_chooser_dialog_new("Save handle configuration", window, GTK_FILE_CHOOSER_ACTION_SAVE, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT, NULL);
        if(gtk_dialog_run(GTK_DIALOG(d)) == GTK_RESPONSE_ACCEPT) {
            const char* filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));
            FILE* f = fopen(filename, "w");
            int l = current_toy->handles.size();
            for(int i = 0; i < l; i++)
                fprintf(f, "%lf %lf\n", current_toy->handles[i][0], current_toy->handles[i][1]);
            fclose(f);
        }
        gtk_widget_destroy(d);
    }
}

void save_cairo() {
    GtkWidget* d = gtk_file_chooser_dialog_new("Save file as svg or pdf", window, GTK_FILE_CHOOSER_ACTION_SAVE, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT, NULL);
    if(gtk_dialog_run(GTK_DIALOG(d)) == GTK_RESPONSE_ACCEPT) {
        const char* filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));
        cairo_surface_t* cr_s;
        int l = strlen(filename);
        #if CAIRO_HAS_PDF_SURFACE
        if (l >= 4 && strcmp(filename + l - 4, ".pdf") == 0)
            cr_s = cairo_pdf_surface_create(filename, 600., 600.);
        #endif
        #if CAIRO_HAS_SVG_SURFACE
        #if CAIRO_HAS_PDF_SURFACE        
        else
        #endif
            cr_s = cairo_svg_surface_create(filename, 600., 600.);
        #endif
        cairo_t* cr = cairo_create(cr_s);
        
        if(current_toy != NULL)
            current_toy->draw(cr, new std::ostringstream, 600, 600, true);

        cairo_show_page(cr);
        cairo_destroy (cr);
        cairo_surface_destroy (cr_s);
    }
    gtk_widget_destroy(d);
}

void save_image() {
    GtkWidget* d = gtk_file_chooser_dialog_new("Save file as png", window, GTK_FILE_CHOOSER_ACTION_SAVE, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT, NULL);
    if(gtk_dialog_run(GTK_DIALOG(d)) == GTK_RESPONSE_ACCEPT) {
        const char* filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));
        cairo_surface_t* cr_s = cairo_image_surface_create ( CAIRO_FORMAT_ARGB32, 600, 600 );
        cairo_t* cr = cairo_create(cr_s);
        
        if(current_toy != NULL)
            current_toy->draw(cr, new std::ostringstream, 600, 600, true);

        cairo_show_page(cr);
        cairo_surface_write_to_png(cr_s, filename);
        cairo_destroy (cr);
        cairo_surface_destroy (cr_s);
    }
    gtk_widget_destroy(d);
}

static gint delete_event(GtkWidget* window, GdkEventAny* e, gpointer data) {
    (void)( window);
    (void)( e);
    (void)( data);

    gtk_main_quit();
    return FALSE;
}

static gboolean expose_event(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
    (void)(data);
    cairo_t *cr = gdk_cairo_create(widget->window);
    
    int width = 256;
    int height = 256;
    gdk_drawable_get_size(widget->window, &width, &height);

    std::ostringstream notify;

    if(current_toy != NULL) current_toy->draw(cr, &notify, width, height, false);
    cairo_destroy(cr);

    return TRUE;
}

static gint mouse_motion_event(GtkWidget* widget, GdkEventMotion* e, gpointer data) {
    (void)(data);
    (void)(widget);

    if(current_toy != NULL) current_toy->mouse_moved(e);

    return FALSE;
}

static gint mouse_event(GtkWidget* widget, GdkEventButton* e, gpointer data) {
    (void)(data);
    (void)(widget);

    if(current_toy != NULL) current_toy->mouse_pressed(e);

    return FALSE;
}

static gint mouse_release_event(GtkWidget* widget, GdkEventButton* e, gpointer data) {
    (void)(data);
    (void)(widget);

    if(current_toy != NULL) current_toy->mouse_released(e);

    return FALSE;
}

static gint key_release_event(GtkWidget *widget, GdkEventKey *e, gpointer data) {
    (void)(data);
    (void)(widget);

    if(current_toy != NULL) current_toy->key_hit(e);

    return FALSE;
}

GtkItemFactoryEntry menu_items[] = {
    { "/_File",             NULL,           NULL,           0,  "<Branch>"                    },
    { "/File/_Open Handles","<CTRL>O",      open,           0,  "<StockItem>", GTK_STOCK_OPEN },
    { "/File/_Save Handles","<CTRL>S",      save,           0,  "<StockItem>", GTK_STOCK_SAVE_AS },
    { "/File/sep",          NULL,           NULL,           0,  "<Separator>"                 },
    { "/File/Save SVG/PDF", NULL,           save_cairo,     0,  "<StockItem>", GTK_STOCK_SAVE },
    { "/File/Save PNG",     NULL,           save_image,     0,  "<StockItem>", GTK_STOCK_SELECT_COLOR }, 
    { "/File/sep",          NULL,           NULL,           0,  "<Separator>"                 },
    { "/File/_Quit",        "<CTRL>Q",      gtk_main_quit,  0,  "<StockItem>", GTK_STOCK_QUIT },
    { "/_Help",             NULL,           NULL,           0,  "<LastBranch>"                },
    { "/Help/About",        NULL,           make_about,     0,  "<StockItem>", GTK_STOCK_ABOUT}
};
gint nmenu_items = 9;

void init(int argc, char **argv, char *title, Toy* t) {
    current_toy = t;
    gtk_init (&argc, &argv);
    
    gdk_rgb_init();

    window = GTK_WINDOW(gtk_window_new(GTK_WINDOW_TOPLEVEL));

    gtk_window_set_title(GTK_WINDOW(window), title);

    //Creates the menu from the menu data above
    GtkAccelGroup* accel_group = gtk_accel_group_new();
    GtkItemFactory* item_factory = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", accel_group);
    gtk_item_factory_create_items(item_factory, nmenu_items, menu_items, NULL);
    gtk_window_add_accel_group(window, accel_group);
    GtkWidget* menu = gtk_item_factory_get_widget(item_factory, "<main>");

    //gtk_window_set_policy(GTK_WINDOW(window), TRUE, TRUE, TRUE);

    gtk_signal_connect(GTK_OBJECT(window), "delete_event", GTK_SIGNAL_FUNC(delete_event), NULL);

    gtk_widget_push_visual(gdk_rgb_get_visual());
    gtk_widget_push_colormap(gdk_rgb_get_cmap());
    canvas = gtk_drawing_area_new();

    gtk_widget_add_events(canvas, (GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | GDK_KEY_PRESS_MASK | GDK_POINTER_MOTION_MASK));

    gtk_signal_connect(GTK_OBJECT (canvas), "expose_event", GTK_SIGNAL_FUNC(expose_event), 0);
    gtk_signal_connect(GTK_OBJECT(canvas), "button_press_event", GTK_SIGNAL_FUNC(mouse_event), 0);
    gtk_signal_connect(GTK_OBJECT (canvas), "button_release_event", GTK_SIGNAL_FUNC(mouse_release_event), 0);
    gtk_signal_connect(GTK_OBJECT (canvas), "motion_notify_event", GTK_SIGNAL_FUNC(mouse_motion_event), 0);
    gtk_signal_connect(GTK_OBJECT(canvas), "key_press_event", GTK_SIGNAL_FUNC(key_release_event), 0);

    gtk_widget_pop_colormap();
    gtk_widget_pop_visual();

    GtkWidget* box = gtk_vbox_new (FALSE, 0);
    gtk_container_add(GTK_CONTAINER(window), box);

    gtk_box_pack_start (GTK_BOX (box), menu, FALSE, FALSE, 0);

    GtkWidget* pain = gtk_vpaned_new();
    gtk_box_pack_start(GTK_BOX(box), pain, TRUE, TRUE, 0);
    gtk_paned_add1(GTK_PANED(pain), canvas);

    gtk_window_set_default_size(GTK_WINDOW(window), 600, 600);

    gtk_widget_show_all(GTK_WIDGET(window));

    // Make sure the canvas can receive key press events.
    GTK_WIDGET_SET_FLAGS(canvas, GTK_CAN_FOCUS);
    assert(GTK_WIDGET_CAN_FOCUS(canvas));
    gtk_widget_grab_focus(canvas);
    assert(gtk_widget_is_focus(canvas));

    gtk_main();
}
