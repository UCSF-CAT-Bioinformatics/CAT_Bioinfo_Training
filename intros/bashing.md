# Course Title:** Bash Scripting in a Flash (with AI Power!)

A 1-hour introductory Bash scripting course, culminating in using AI assistance:

**Target Audience:** Beginners with basic computer literacy.

**Course Objectives:**

*   Understand what Bash scripting is and why it's useful.
*   Learn fundamental Bash commands and syntax.
*   Write simple Bash scripts.
*   Discover how AI tools like Gemini can aid in script development.


## Introduction & What is Bash? (5 minutes) 

A command-line interpreter and scripting language for Unix-like systems (Linux, macOS).

Highlight its uses: automating tasks, system administration, file manipulation, etc.

Knowing how to use Bash and script in bash provides give you power over your system and continues to be relevant even in today's GUI-driven world.

## Recalling Basic Bash Commands & Syntax: (15 minutes) 

*   **Navigation:**
    `pwd` (print working directory), `cd` (change directory), `ls` (list files).
*   **File Manipulation:**
    `touch` (create file), `mkdir` (make directory), `rm` (remove), `cp` (copy), `mv` (move).
*   **Viewing Files:**
*   `cat` (concatenate and display), `less` (page through), `head` (first few lines), `tail` (last few lines).
*   **Input/Output Redirection:**
*   `>` (output to file, overwriting), `>>` (append to file), `<` (input from file), `|` (pipe output to another command).
*   **Example:**
*   Create a directory, create a file inside it, write some text to the file (echo, nano), and display its contents.

## Writing a Simple Bash Script: (20 minutes) 

*   **Shebang:**
*   Explain the `#!/bin/bash` line at the beginning of a script (tells the system which interpreter to use).
*   **Variables:**
*   Demonstrate how to declare and use variables (e.g., `name="John"`, `echo "Hello, $name"`).
*   **Echo:**
*   Explain the `echo` command for displaying output.
*   **Comments:**
*   Show how to add comments using `#`.
*   **Permissions:**
*   Explain `chmod +x script.sh` to make the script executable.

### **Simple Script Example:**

```bash
#!/bin/bash
name="User"
echo "Hello, $name!"
echo "Today is $(date)"
```
*   **Walk through the script:** Explain each line and its purpose.
*   **Hands-on exercise:** Have participants create a simple script that prints their name and the current date.

## Using AI (Gemini) for Bash Scripting (15 minutes) 

*   **Introduce Gemini:** Explain how AI models like Gemini can assist in coding.
*   **Demonstrate prompting Gemini:**
    *   **Simple Task:** "Write a bash script to list all files in the current directory that are larger than 1MB."
    *   **More Complex Task:** "Write a bash script that takes a directory as input and counts the number of files of each type (e.g., .txt, .pdf, .jpg) in that directory."
*   **Analyze Gemini's output:**
    *   Discuss the generated code, explaining its logic and syntax.
    *   Show how to test and modify the code.
    *   Emphasize the importance of understanding the generated code, not just blindly copying it.
*   **Prompt Engineering Tips for Gemini:**
    *   Be specific and descriptive in your prompts.
    *   Provide context and examples if needed.
    *   Iterate on your prompts based on the results.
*   **Benefits of using AI for scripting:** Speed, discovering new methods, learning best practices.
*   **Limitations:** AI may not always produce perfect code and requires human review.

**Materials:**

*   Slides (optional, but recommended for visual aids).
*   Text editor (e.g., Nano, Vim, VS Code) for demonstrations.
*   Access to a terminal (Linux/macOS or Windows Subsystem for Linux (WSL)).
*   Internet access for demonstrating Gemini.
