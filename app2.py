app.pyfrom shiny import App, render, ui

# 1. Define the User Interface (UI)
app_ui = ui.page_fluid(
    ui.h2("Simple Shiny Core App"),
    ui.input_text("name", "Enter your name:", value="World"),
    ui.output_text_verso("greeting") # or output_text
)

# Alternative standard output component for text
app_ui = ui.page_fluid(
    ui.h2("Simple Shiny Core App"),
    ui.input_text("name", "Enter your name:", value="World"),
    ui.output_text("greeting")
)

# 2. Define the Server Logic
def server(input, output, session):
    @output
    @render.text
    def greeting():
        return f"Hello, {input.name()}!"

# 3. Combine UI and Server into an App instance
app = App(app_ui, server)
