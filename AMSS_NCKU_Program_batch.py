
##################################################################
##
## AMSS-NCKU 无交互（批处理）启动程序 | AMSS-NCKU non-interactive (batch) startup program
## 用于超算作业脚本等无法手动输入的环境 | Intended for environments such as supercomputer job scripts where manual input is impossible
##
##################################################################

import builtins
import os
import runpy
import sys


def select_language():
    """从环境变量 AMSS_NCKU_LANGUAGE 读取输出语言，默认 Chinese。
    Reads the output language from the environment variable AMSS_NCKU_LANGUAGE; defaults to Chinese.

    可选值: Chinese、English（大小写不敏感）。
    Allowed values: Chinese, English (case-insensitive).
    """
    language = os.environ.get("AMSS_NCKU_LANGUAGE", "Chinese").strip().lower()
    if language in ("english", "en"):
        return "English"
    return "Chinese"


def continue_or_stop( File_directory, input_language ):
    """无交互版本：输出目录存在时自动退出，不再询问用户，也不覆盖。 | Non-interactive version: exits automatically if the output directory exists, without asking the user again or overwriting it."""
    if os.path.exists(File_directory):
        if ( input_language == "Chinese" ):
            print(                                                                               )
            print( " 无交互模式：输出目录已存在，为避免覆盖自动退出计算 ", File_directory          )
            print( " 如需重新计算，请先删除该目录，或在输入文件中更换输出目录名称 "                 )
            print(                                                                               )
        elif ( input_language == "English" ):
            print(                                                                                                    )
            print( " Non-interactive mode: output directory exists, exiting automatically to avoid overwriting ", File_directory )
            print( " To rerun the simulation, remove the directory or change the output directory in the input file "          )
            print(                                                                                                    )
        sys.exit( 1 )
    return


def main():

    program_directory = os.path.dirname(os.path.abspath(__file__))
    if program_directory not in sys.path:
        sys.path.insert(0, program_directory)

    import main_program_function

    input_language = select_language()

    if ( input_language == "Chinese" ):
        print(                                                                                       )
        print( " 无交互模式：屏幕输出语言 =", input_language, "（可用环境变量 AMSS_NCKU_LANGUAGE 修改）" )
        print(                                                                                       )
    elif ( input_language == "English" ):
        print(                                                                                                  )
        print( " Non-interactive mode: language =", input_language, "（set AMSS_NCKU_LANGUAGE to change）"       )
        print(                                                                                                  )

    ## 固定语言，跳过交互式语言选择 | Fix the language, skipping the interactive language selection
    main_program_function.setting_language = lambda: input_language

    ## 输出目录存在时自动退出，跳过交互式询问 | Exit automatically when the output directory exists, skipping the interactive prompt
    main_program_function.continue_or_stop = continue_or_stop

    ## 所有"按回车继续"的提示都直接跳过 | All "press Enter to continue" prompts are skipped directly
    builtins.input = lambda *args, **kwargs: ""

    runpy.run_path( os.path.join(program_directory, "AMSS_NCKU_Program.py"), run_name="__main__" )


if __name__ == "__main__":
    main()
